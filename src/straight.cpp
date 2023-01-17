//
// Created by christoph on 14.01.23.
//

#include "straight.h"
#include <iostream>

/// private
void straight_generator::calculate_mass_center() {
    // add all the points up, divide by the number of points
    mass_center.setZero();
    for(auto point : points->boundary_points) {
        mass_center += point->point;
    }
    mass_center /= double(points->boundary_points.size());
}

void straight_generator::calculate_all_straights() {
    auto iter = points->boundary_points.begin();
    while(iter != points->boundary_points.end()) {
        // special cases for end, we set the next iter to iter + 1
        auto next_iter = points->boundary_points.begin();
        // check if not the last one
        if(!((iter+1) == points->boundary_points.end())) {
            next_iter = iter + 1;
        }
        // fill the values
        auto s  = new surface_t ;
        s->point = (next_iter.operator*()->point - iter.operator*()->point);
        // rotate the vector by 90 degrees forward (doesnt really matter which direction)
        s->direction = {s->point.y(),-s->point.x()};
        // put the point in the middle
        s->point /= 0.5;
        surfaces.push_back(s);
        iter++;
    }
}


/// public
straight_generator::straight_generator(boundaryPointConstructor *p) {
    points = p;
}

void straight_generator::init() {
    calculate_mass_center();
    calculate_all_straights();
}

int straight_generator::calculate_intersections(nodePoint_t* point) {
    int number_of_intersections = 0;
    // determine straight to the mass center


    return number_of_intersections;
}
