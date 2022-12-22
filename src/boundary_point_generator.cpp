
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
    //
    boundaryType_t type = BOUNCE_BACK;
    point_t current;
    current.setZero();
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