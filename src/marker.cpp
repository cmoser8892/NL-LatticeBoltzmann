#include "marker.h"
#include <fstream>
#include <iostream>
/**
 * Constructor with lots of optional fields.
 * @param s
 * @param kd
 * @param md
 */
markerIBM::markerIBM(straightGenerator *s, double md) {
    sg = s;
    marker_distance = md;
}

/**
 * De-constructor gets rid of the points.
 */
markerIBM::~markerIBM() {
    // delete all the info in the vector
    for(auto p : marker_points) {
        delete p;
    }
}

/**
 * Marker distributor, main goal is too evenly distribute the markers on the surface.
 * @attention surface is blind of its own inner structure
 */
void markerIBM::distribute_markers() {
    // start at the first surface and look
    // todo do i need to know which is the direction towards the fluid?!
    if(sg != nullptr) {
        double total_surface = sg->calculate_total_surface_length();
        // the markers should be equally space over the total surface
        // we dont want any rest we just want equal spacing
        // determine the number of markers in a surface
        double markers_fit_in= std::floor(total_surface/marker_distance);
        // check weather or not we have to equalize so that everything is equally spaced
        // watch out for low marker numbers
        double add_on = std::fmod(total_surface,marker_distance);
        marker_distance += add_on/markers_fit_in;
        // we now step through the surface again
        // setup walk distance initial
        double walk_distance = marker_distance;
        for(auto s : sg->surfaces) {
            // in all directions just walk the distance
            vector_t walker = s->direction.normalized() * walk_distance;
            point_t marker = s->point;
            bool end = false;
            double overshoot = 0;
            // go through the individual surface
            while(!end) {
                walker = s->direction.normalized() * walk_distance;
                marker += walker;
                // still on
                if(point_on_straight(s,&marker,&overshoot)) {
                    // add
                    walk_distance = marker_distance;
                    // put the point into a permanent object
                    auto add_me = new point_t;
                    *add_me = marker;
                    marker_points.push_back(add_me);
                }
                // overshoot
                else {
                    // overshoot processing
                    walk_distance = overshoot;
                    end = true;
                }
            }
        }
    }
}

void markerIBM::individual_distribute_markers(straight_t s) {
    double length = s.max_t;
    // determine the number of markers in a surface
    double markers_fit_in= std::floor(length/marker_distance);
    // check weather or not we have to equalize so that everything is equally spaced
    // watch out for low marker numbers
    double add_on = std::fmod(length,marker_distance);
    marker_distance += add_on/markers_fit_in;
    double walk_distance = marker_distance;
    // in all directions just walk the distance
    vector_t walker = s.direction.normalized() * walk_distance;
    point_t marker = s.point;
    bool end = false;
    double dk = 0;
    // go through the individual surface
    while(!end) {
        walker = s.direction.normalized() * walk_distance;
        // still on
        if(point_on_straight(&s,&marker,&dk)) {
            // add
            walk_distance = marker_distance;
            // put the point into a permanent object
            auto add_me = new point_t;
            *add_me = marker;
            marker_points.push_back(add_me);
        }
        // overshoot
        else {
            // overshoot processing
            end = true;
        }
        marker += walker;
    }
}

/**
 * Gives out the actual marker distance used.
 * @return
 */
double markerIBM::return_marker_distance() {
    return marker_distance;
}

/**
 * Writes out the marker points into a file or not.
 * @param write_file
 */
void markerIBM::write_out_markers(bool write_file) {
    std::ofstream out;
    out.open("markers");
    if(out.is_open()) {
        for(auto m : marker_points) {
            out << m->x() << " " << m->y() << std::endl;
        }
    }
}

void markerIBM::set_marker_distance(double md) {
    marker_distance = md;
}