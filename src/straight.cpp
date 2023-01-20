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
        auto s  = new surface_t;
        s->point = iter.operator*()->point;
        // rotate the vector by 90 degrees forward (doesnt really matter which direction)
        s->direction = (next_iter.operator*()->point - iter.operator*()->point);
        // put the point in the middle
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

int straight_generator::calculate_intersections(nodePoint_t* node_point) {
    /**
     * 3 passes have to be made to calculate to calcuate a valid intersection
     *  1 does the straight hit the surface in the area between the two points that define it
     *  2 how does the straight hit the surface (posetive or negative we only care about posetiv
     *  3 have we already hit an edgepoint
     */
    int number_of_intersections = 0;
    // corner case mass center lays on point
    point_t p = node_point->position;
    if(p == mass_center) {
        return 1;
    }
    // determine straight to the mass center
    straight_t straight;
    straight.point = node_point->position; // => r
    straight.direction = mass_center - straight.point;
    vector_t normal = {straight.direction.y(), -straight.direction.x()}; // => n
    // go through the surface and take a look
    std::vector<double> previous_ss;
    for(auto surf : surfaces) {
        // t = ((r - o)·n)/(n·d)
        // surf->point => o
        // surf->direction => d
        double t = ((straight.point - surf->point).dot(normal))/
                   (normal.dot(surf->direction));
        if((t >= 0.0) && (t <= 1.0)) {
            // check if direction of the finding is posetiv in the direction of the vector
            vector_t surface_normal = {surf->direction.y(),-surf->direction.x()};
            double s = ((surf->point-straight.point).dot(surface_normal))/(surface_normal.dot(straight.direction));
            if(s >= 0) {
                bool add = true;
                // might be a bad idea no just check s and not o +sd
                // but it should be alright
                for (auto ps : previous_ss) {
                    if(ps == s) {
                        add = false;
                    }
                }
                previous_ss.push_back(s);
                if(add) {
                    number_of_intersections++;
                }
            }
        }
    }
    return number_of_intersections;
}

bool straight_generator::node_inside(nodePoint_t *point) {
    // even out; odd in
    int value = calculate_intersections(point);
    return ((value%2) == 0);
}
