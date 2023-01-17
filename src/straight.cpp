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
        auto s  = new straight_t;
        s->point = iter.operator*()->point;
        s->direction = (next_iter.operator*()->point - iter.operator*()->point);
        s->length = s->direction.norm();
        straights.push_back(s);
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

void straight_generator::calculate_intersections() {
    // todo
}
