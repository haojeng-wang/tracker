#ifndef __TRACKER_H__
#define __TRACKER_H__

#include <stdint.h>
#include <deque>
#include <vector>

#include "scv.h"

using namespace scv;
using namespace std;

#define PARTICLE_NUMBER 100 
#define HIST_BIN_NUM    4096

class particle
{

public:
    particle() {}

    particle(
        int _x, 
        int _y, 
        int _prex, 
        int _prey, 
        int _rect_x,
        int _rect_y,
        int _rect_width,
        int _rect_height,
        float _weight
        ) : 
        x(_x), 
        y(_y), 
        prex(_prex), 
        prey(_prey), 
        rect_x(_rect_x),
        rect_y(_rect_y),
        rect_width(_rect_width),
        rect_height(_rect_height),
        weight(_weight)
    {

    }

    //particle(const particle& p)
    //{
    //    x = p.x;
    //    y = p.y;
    //    prex = p.prex;
    //    prey = p.prey;
    //    rect_x = p.rect_x;
    //    rect_y = p.rect_y;
    //    rect_width = p.rect_width;
    //    rect_height = p.rect_height;
    //    weight = p.weight;
    //}

    int x, y;
    int prex, prey;
    int rect_x, rect_y, rect_width, rect_height;
    float weight;

};


class tracker
{
public:
    tracker() {}

    tracker(int sel_x, int sel_y, int sel_w, int sel_h, bool sel_acti = true) : tr_x(sel_x), tr_y(sel_y), tr_w(sel_w), tr_h(sel_h), activated(sel_acti)
    {           
    }

    tracker(const tracker & tra) 
    {
        tr_x = tra.tr_x;
        tr_y = tra.tr_y;
        tr_w = tra.tr_w;
        tr_h = tra.tr_h;        
    }

    ~tracker() {}

    virtual void init(const Matrix& img) = 0;
    virtual void track(const Matrix& img) = 0;

    Rect get_rect() const 
    {
        return Rect(tr_x, tr_y, tr_w, tr_h);
    }

protected:
    uint32_t tr_x;
    uint32_t tr_y;
    uint32_t tr_w;
    uint32_t tr_h;

    bool activated;

};

class pftracker : public tracker
{
public:
    pftracker() {};

    pftracker(int sel_x, int sel_y, int sel_w, int sel_h) : tracker(sel_x, sel_y, sel_w, sel_h, true)
    {     
        particles.resize(PARTICLE_NUMBER, particle(
            sel_x + (sel_w >> 1),
            sel_y + (sel_h >> 1),
            sel_x + (sel_w >> 1),
            sel_y + (sel_h >> 1),
            sel_x,
            sel_y,
            sel_w,
            sel_h,
            0.0f));        
    }

    //Copy constructor
    pftracker(const pftracker& pftra) : tracker(pftra)
    {
        memcpy(hist, pftra.hist, HIST_BIN_NUM*sizeof(uint32_t));
        particles = pftra.particles;        
    }

    ~pftracker()
    {
        particles.clear();        
    }

    void init(const Matrix& img);
    void track(const Matrix& img);

private:
    uint32_t hist [HIST_BIN_NUM];
    vector<particle> particles;

};




#endif 