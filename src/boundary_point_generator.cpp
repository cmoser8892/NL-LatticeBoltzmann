#include <iostream>
#include "boundary_point_generator.h"

boundaryPointConstructor::boundaryPointConstructor(point_t s) {
    size = s;
    limits << s.x()-1,s.y()-1;
}

void boundaryPointConstructor::one_direction(int limit, vector_t dir, point_t *start, boundaryType_t  b) {
    // goes up and down a line and inizlizes it
    for(int i = 0; i < limit; i++) {
        set_point(start,b);
        *start += dir;
    }
}

void boundaryPointConstructor::set_point(point_t* p, boundaryType_t b) {
    auto boundary_point = new boundaryPoint_t;
    boundary_point->point = *p;
    boundary_point->type = b;
    boundary_points.push_back(boundary_point);
}

void boundaryPointConstructor::init_quader() {
    point_t current;
    current.setZero();
    init_quader(current);
}

void boundaryPointConstructor::init_chopped_quader(point_t point, int chopfactor) {
    if(chopfactor == 0) {
        chopfactor = 2147483647;
    }
    if(chopfactor < 2) {
        throw std::runtime_error("unrealitic chopfactor");
    }
    boundaryType_t type = BOUNCE_BACK;
    point_t current = point;
    // go from 0 till the in x directions
    int chop_x = int(size.x()/chopfactor);
    int chop_y = int(size.y()/chopfactor);
    vector_t direction;
    // go through x
    direction = {1,0};
    one_direction(int(limits.x())-chop_x,direction,&current, type);
    direction = {0,1};
    one_direction(chop_y,direction,&current, type);
    direction = {1,0};
    one_direction(chop_x,direction,&current, type);
    direction = {0,1};
    one_direction(int(limits.y())-chop_y,direction,&current, type);
    // go through y
    direction = {-1,0};
    one_direction(int(limits.x()),direction,&current, type);
    // go through x
    direction = {0,-1};
    one_direction(int(limits.y()),direction,&current, type);
}

void boundaryPointConstructor::init_quader(point_t p) {
    //
    boundaryType_t type = BOUNCE_BACK;
    point_t current = p;
    // go from 0 till the in x directions
    vector_t direction;
    // go through x
    direction = {1,0};
    one_direction(int(limits.x()),direction,&current, type);
    // go through x
    direction = {0,1};
    one_direction(int(limits.y()),direction,&current, type);
    // go through x
    direction = {-1,0};
    one_direction(int(limits.x()),direction,&current, type);
    // go through x
    direction = {0,-1};
    one_direction(int(limits.y()),direction,&current, type);
}

void boundaryPointConstructor::init_sliding_lid() {
    // greate a slinding lid container with the given sizes
    //in our case y max is the boundary that is moving
    // so all boundaries with y = y_max are BOUNDARY_MOVING
    // init a quader
    double limit_y = limits.y();
    init_quader();
    for(auto b : boundary_points) {
        if(b->point.y() == limit_y) {
            b->type = BOUNCE_BACK_MOVING;
        }
    }
}


void boundaryPointConstructor::init_chopped_sliding_lid(point_t start, int chopfactor) {
    double limit_y = limits.y() + start.y();
    init_chopped_quader(start,chopfactor);
    for(auto b : boundary_points) {
        if(b->point.y() == limit_y) {
            b->type = BOUNCE_BACK_MOVING;
        }
    }
}