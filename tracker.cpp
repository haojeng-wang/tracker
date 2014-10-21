#include <random>
#include <cmath>

#include "tracker.h"
#include "scv.h"


const int A1 = 2;
const int A2 = -1;
const int B0 = 1;
const float sigmax = 8.0f;
const float sigmay = 8.0f;
const uint32_t standard_size = 5; //2^x = 32 (ex: 2^5 = 32)
std::default_random_engine generator;

inline bool struct_cmp_by_freq(const particle& a, const particle& b)
{
    return a.weight > b.weight;
}

inline void generateHist(const Matrix& select_image, uint32_t* hist_range, const uint32_t& hist_size, const Rect& rt)
{ 
    memset((void*)hist_range, 0, hist_size*sizeof(uint32_t));

    int index_base, index_base_i;
    int h_ind, s_ind, v_ind, hsv_ind;

    uint8_t* data_hsv = (uint8_t*)select_image.data;

    int x = rt.x;
    int y = rt.y;
    int final_x = x + rt.width;
    int final_y = y + rt.height; 

    int width_step = 3*select_image.cols;

    for(int i = y; i < final_y; ++i) {

        index_base_i = i*width_step;

        for(int j = x; j < final_x; ++j) {

            index_base = j*3 + index_base_i;

            //h_ind = min(15, int(data_hsv[index_base + 0]/11.25f));
            h_ind = min(15, data_hsv[index_base + 0] >> 4);  //16            
            s_ind = min(15, data_hsv[index_base + 1] >> 4);  //16     
            v_ind = min(15, data_hsv[index_base + 2] >> 4);  //16      

            hsv_ind = (h_ind << 8) + (s_ind << 4) + v_ind;

            hist_range[hsv_ind] += 1;

        }
    }
}

inline float compareHistBhattacharyya(const uint32_t* hist1, const uint32_t* hist2, const uint32_t& hist_size)
{
    float sum_hist1 = 0, sum_hist2 = 0, sum_product_sqrt = 0;
    float value1, value2;
    for(uint32_t i = 0; i < hist_size; ++i) {
        value1 = (float)hist1[i];
        value2 = (float)hist2[i];
        sum_hist1 += value1;
        sum_hist2 += value2;
        sum_product_sqrt += sqrt(value1 * value2);
    }
    return sqrt(1 - sum_product_sqrt/sqrt(sum_hist1 * sum_hist2));
    
}


void pftracker::init(const Matrix& img)
{   
    generateHist(img, hist, HIST_BIN_NUM, Rect(tr_x, tr_y, tr_w, tr_h));
}

void pftracker::track(const Matrix& img)
{
    const uint32_t cam_w = img.cols;
    const uint32_t cam_h = img.rows;

    /****Upadte particle****/
    particle* p_particle = NULL;

    float sum = 0.0f;

    const int sigmax_gain = (tr_w >> standard_size);
    const int sigmay_gain = (tr_h >> standard_size);
    const float sigx = sigmax * sigmax_gain;
    const float sigy = sigmay * sigmay_gain;

    std::normal_distribution<float> distrix(0, sigx);
    std::normal_distribution<float> distriy(0, sigy);

    uint32_t track_hist [HIST_BIN_NUM];

    /****Upadte particle****/
    for(int i = 0; i < PARTICLE_NUMBER; ++i)
    {
        int x, y;
        int xpre, ypre;

        p_particle = &particles[i];

        xpre = p_particle->x;
        ypre = p_particle->y;

        float dx = distrix(generator); //Generate gaussian random variable
        float dy = distriy(generator); //Generate gaussian random variable

        x = (int)(A1 * (p_particle->x) + A2 * (p_particle->prex) + B0 * dx);
        y = (int)(A1 * (p_particle->y) + A2 * (p_particle->prey) + B0 * dy);

        p_particle->x = max(0, min(x, int32_t(cam_w - 1)));
        p_particle->y = max(0, min(y, int32_t(cam_h - 1)));

        p_particle->prex = xpre;
        p_particle->prey = ypre;

        p_particle->rect_x	    = max(0, min((p_particle->x - (p_particle->rect_width  >> 1)), int32_t(cam_w) - 1));
        p_particle->rect_y      = max(0, min((p_particle->y - (p_particle->rect_height >> 1)), int32_t(cam_h) - 1));
        p_particle->rect_width  = min((p_particle->rect_width ), int32_t(cam_w) - 1 - p_particle->rect_x);
        p_particle->rect_height = min((p_particle->rect_height), int32_t(cam_h) - 1 - p_particle->rect_y);

        generateHist(img, track_hist, HIST_BIN_NUM, Rect(p_particle->rect_x, p_particle->rect_y, p_particle->rect_width, p_particle->rect_height));	

        p_particle->weight = 1.0f - compareHistBhattacharyya(hist, track_hist, HIST_BIN_NUM);

        sum += p_particle->weight;

    }

    /****Normalization****/
    for(int i = 0; i < PARTICLE_NUMBER; ++i)
        particles[i].weight /= sum;

    sort(particles.begin(), particles.end(), struct_cmp_by_freq);

    /****Resample****/    
    vector<particle> newparticles(PARTICLE_NUMBER);
    particle* p_new_particle = NULL;
    int np = 0, k = 0;

    for(int i = 0; i < PARTICLE_NUMBER; ++i)
    {
        p_particle = &particles[i];

        np = (int)(p_particle->weight*PARTICLE_NUMBER);

        for(int j = 0; j < np; ++j)
        {
            p_new_particle = &newparticles[k];

            p_new_particle->x          = p_particle->x;
            p_new_particle->y          = p_particle->y;
            p_new_particle->prex	   = p_particle->prex;
            p_new_particle->prey	   = p_particle->prey;
            p_new_particle->rect_x     = p_particle->rect_x;
            p_new_particle->rect_y     = p_particle->rect_y;
            p_new_particle->rect_width = p_particle->rect_width;
            p_new_particle->rect_height= p_particle->rect_height;
            p_new_particle->weight     = p_particle->weight;

            k++;

            if( k == PARTICLE_NUMBER)
                goto EXITOUT;
        }

    }       

    p_particle = &particles[0];

    while( k < PARTICLE_NUMBER)
    {        
        p_new_particle = &newparticles[k];

        p_new_particle->x           = p_particle->x;
        p_new_particle->y           = p_particle->y;
        p_new_particle->prex	    = p_particle->prex;
        p_new_particle->prey	    = p_particle->prey;
        p_new_particle->rect_x      = p_particle->rect_x;
        p_new_particle->rect_y      = p_particle->rect_y;
        p_new_particle->rect_width  = p_particle->rect_width;
        p_new_particle->rect_height = p_particle->rect_height;
        p_new_particle->weight      = p_particle->weight;

        k++;
    }

EXITOUT:
    for(int i = 0; i < PARTICLE_NUMBER; ++i)
    {
        p_particle     = &particles[i];
        p_new_particle = &newparticles[i];

        p_particle->x          = p_new_particle->x;
        p_particle->y          = p_new_particle->y;
        p_particle->prex       = p_new_particle->prex;
        p_particle->prey       = p_new_particle->prey;
        p_particle->rect_x     = p_new_particle->rect_x;
        p_particle->rect_y     = p_new_particle->rect_y;
        p_particle->rect_width = p_new_particle->rect_width;
        p_particle->rect_height= p_new_particle->rect_height;
        p_particle->weight     = p_new_particle->weight;
    }

    float weight_sum = 0.0f;
    float weight_ratio = 0.0f;

    for(int i = 0; i < PARTICLE_NUMBER; ++i)
        weight_sum += particles[i].weight;

    if(weight_sum != 0) {

        float rectTracking_x = 0.0;
        float rectTracking_y = 0.0;
        float rectTracking_w = 0.0;
        float rectTracking_h = 0.0;

        for(int i = 0; i < PARTICLE_NUMBER; ++i)
        {
            p_particle = &particles[i];

            weight_ratio = p_particle->weight/weight_sum;

            rectTracking_x += p_particle->rect_x      * weight_ratio;
            rectTracking_y += p_particle->rect_y      * weight_ratio;
            rectTracking_w += p_particle->rect_width  * weight_ratio;
            rectTracking_h += p_particle->rect_height * weight_ratio;
        }

        tr_x = (uint32_t)rectTracking_x;
        tr_y = (uint32_t)rectTracking_y;
        tr_w = (uint32_t)rectTracking_w;
        tr_h = (uint32_t)rectTracking_h;

    }
}

