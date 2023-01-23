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
    /// surface based algorithm to calcualte intersections
    /**
     * 3 passes have to be made to calculate to calcuate a valid intersection
     *  1 does the straight hit the surface in the area between the two points that define it
     *  2 how does the straight hit the surface (posetive or negative we only care about posetiv
     *  3 have we already hit an edgepoint
     */
    int number_of_intersections = 0;
    // corner case mass center lays on point
    point_t p = node_point->position;
    // determine straight to the mass center
    straight_t straight;
    straight.point = node_point->position; // => r
    straight.direction =  mass_center - straight.point;
    vector_t normal = {straight.direction.y(), -straight.direction.x()}; // => n
    // go through the surface and take a look
    std::vector<point_t> already_found;
    already_found.clear();
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
                point_t point = straight.point + s*straight.direction;
                bool add = true;
                // might be a bad idea no just check s and not o +sd
                // but it should be alright
                for (auto ps : already_found) {
                    if(ps == point) {
                        add = false;
                    }
                }
                already_found.push_back(point);
                if(add) {
                    number_of_intersections++;
                }
            }
        }
    }
    if(0) {
        std::cout << "Result" << std::endl;
        std::cout << node_point->position.x() << " " << node_point->position.y() << std::endl;
        std::cout << number_of_intersections << std::endl;
    }
    return number_of_intersections;
}

bool straight_generator::node_inside(nodePoint_t *point) {
    // even out; odd in
    /// uses a surface representation to calculate weather nodes are inside or outside
    int value = calculate_intersections(point);
    return ((value%2) == 0);
}

bool straight_generator::node_intersections(nodePoint_t *node_point) {
    // doesnt work?!?!?!
    /// written after
    // https://en.wikipedia.org/wiki/Line%E2%80%93line_intersection
    int number_of_intersections = 0;
    double x1 = node_point->position.x();
    double y1 = node_point->position.y();
    double x2 = mass_center.x();
    double y2 = mass_center.y();
    std::vector<point_t> already_found;
    for(auto surf : surfaces) {
        double x3 = surf->point.x();
        double y3 = surf->point.y();
        double x4 = surf->point.x() + surf->direction.x();
        double y4 = surf->point.y() + surf->direction.y();
        // calculate t and u
        double denominator = (x1 - x2) * (y3 - y4) - (y1 - y2) * (x3 - x4);
        double numerator = (x1 - x3) * (y3 - y4) - (y1 - y3) * (x3 - x4);
        //
        double t = numerator / denominator;
        double u = -((x1 - x2) * (y1 - y3) - (y1 - y2) * (x1 - x3)) / denominator;
        if ((1.0 >= t >= 0.0) && (u >= 0)) {
            double x = x1 + t * (x2 - x1);
            double y = y1 + t * (y2 - y1);
            point_t point = {x,y};
            bool add = true;
            for (auto ps : already_found) {
                if(ps == point) {
                    add = false;
                }
            }
            already_found.push_back(point);
            if(add) {
                number_of_intersections++;
            }
        }

    }
    if(1) {
        std::cout << "Result" << std::endl;
        std::cout << node_point->position.x() << " " << node_point->position.y() << std::endl;
        std::cout << number_of_intersections << std::endl;
    }
    return ((number_of_intersections%2) == 0);
}
