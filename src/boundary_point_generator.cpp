
#include "boundary_point_generator.h"

boundaryPointConstructor::boundaryPointConstructor(point_t s) {
    size = s;
    limits << s.x()-1,s.y()-1;
}

void boundaryPointConstructor::one_direction(int limit, vector_t dir, point_t *start, nodeIdentifier_t n) {
    // goes up and down a line and inizlizes it
    for(int i = 0; i < limit; i++) {
        auto boundary_point = new boundaryPoint_t;
        boundary_point->point = *start;
        boundary_point->type = n;
        boundary_points.push_back(boundary_point);
        *start += dir;
    }
}

void boundaryPointConstructor::init_quader() {
    //
    nodeIdentifier_t type = BOUNDARY;
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