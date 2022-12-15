#include "nodeGenerator.h"
#include "helper_functions.h"
#include "simulation.h"
#include <iostream>

nodeGenerator::nodeGenerator(boundaryPointConstructor *p) {
    points = p;
}

void nodeGenerator::linear_generation() {
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

bool nodeGenerator::check_other_boundary_hit(boundaryPoint_t* p,point_t &check_point){
    // returns true if hit or outside false if not
    bool return_value = false;
    point_t check = check_point.base();
    // right now just does a linear search to check if point hit or not
    for(auto point : points->boundary_points) {
        if(compare_two_points(&check,&point->point)) {
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

void nodeGenerator::determine_neighbors() {
    // go in all directions and search for a match then put channel and and respective handle down
    for(auto n : node_infos) {
        // go though relevant channels
        for(int i = 1; i < CHANNELS; ++i) {
            point_t current = n->position + velocity_set.col(i);
            for(auto search : node_infos) {
                auto temp = point_t(search->position);
                // compare with all the other nodes
                if(compare_two_points(&current, &temp)) {
                    // always include all neighbours if we are a wet node
                    bool add_me = false;
                    if(n->type == WET) {
                        add_me = true;
                    }
                    // if we are a dry node only include nodes that are wet
                    else if(n->type == DRY) {
                        if(search->type == WET) {
                            add_me = true;
                        }
                    }
                    // if one of the above conditions holds add
                    if(add_me) {
                        auto link = new toLinks_t;
                        link->handle = search->handle;
                        link->channel = i;
                        n->links.push_back(link);
                    }
                }
            }
        }
    }
}



bool nodeGenerator::read_data_from_file() {
    return false;
}

void nodeGenerator::write_data_to_file() {
    // nop
}

// public
void nodeGenerator::init() {
    if(!read_data_from_file()) {
        linear_generation();
        determine_neighbors();
    }
}