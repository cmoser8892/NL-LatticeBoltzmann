#include "marker.h"

markerIBM::markerIBM(straightGenerator *s, int kd, double md) {
    sg = s;
    kernel_length = kd;
    marker_distance = md;
}

markerIBM::~markerIBM() {
    // delete all the info in the vector
    for(auto p : marker_points) {
        delete p;
    }
}

void markerIBM::distribute_markers() {
    // start at the first surface and look
    // todo do i need to know which is the direction towards the fluid?!
    if(sg != nullptr) {
        double total_surface = 0;
        for(auto s : sg->surfaces) {
            total_surface += s->direction.norm() * s->max_t;
        }
        // the markers should be equally space over the total surface
        // we dont want any rest we just want equal spacing
        // todo prob bad that surface is structure blind
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

void markerIBM::find_neighborhood() {

}

double markerIBM::return_marker_distance() {
    return marker_distance;
}
