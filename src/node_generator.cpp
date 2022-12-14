#include "node_generator.h"
#include "helper_functions.h"
#include <iostream>

node_generator::node_generator(boundaryPointConstructor *p) {
    points = p;
}

void node_generator::linear_generation() {
    int handle_counter = 1;
    // go throu the boundary points starting at a b point and go throu while still discovering new ones
    vector_t one_zero = {1,0}; // is pretty much arbitray
    for(auto p : points->boundary_points) {
        point_t current = p->point + one_zero;
        while(!check_other_boundary_hit(p,current)) {
            auto n = new nodePoint_t;
            n->handle = handle_counter;
            n->position = current;
            n->type = WET;
            n->boundary = NO_BOUNDARY;
            // dont forget to increase the handle counter each time
            handle_counter++;
            node_infos.push_back(n);
            current += one_zero;
        }
    }
    // lastlly add the boundary points
    for(auto p : points->boundary_points) {
        // std::cout << p << std::endl;
        // set all the variables except the links
        auto n = new nodePoint_t;
        n->handle = handle_counter;
        n->position = p->point;
        n->type = DRY;
        n->boundary = p->type;
        // dont forget to increase the handle counter each time
        handle_counter++;
        node_infos.push_back(n);
    }
}

bool node_generator::check_other_boundary_hit(boundaryPoint_t* p,point_t &check_point){
    // returns true if hit or outside false if not
    bool return_value = false;
    point_t check = check_point.base();
    // right now just does a linear search to check if point hit or not
    for(auto p : points->boundary_points) {
        if(compare_two_points(&check,&p->point)) {
            return_value = true;
            break;
        }
    }
    // then check if still inside sim space
    point_t limit_lower = {0,0};
    point_t limit_upper = points->limits;
    if(!check_inside_limits_upper_lower(&check,&limit_lower,&limit_upper)) {
        return_value = true;
    }
    return return_value;
}

void node_generator::determine_neighbors() {
    // nop
}

void node_generator::init() {
    linear_generation();
    determine_neighbors();
}