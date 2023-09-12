#include "marker.h"
#include <fstream>
#include <iostream>

// private
void markerPoints::walking(straight_t *s, double walker_dist) {

}

// public
/**
 * Constructor with lots of optional fields.
 * @param s
 * @param kd
 * @param md
 */
markerPoints::markerPoints(straightGenerator *s, double md) {
    sg = s;
    marker_distance = md;
}

/**
 * De-constructor gets rid of the points.
 */
markerPoints::~markerPoints() {
    // delete all the info in the vector
    for(auto p : marker_points) {
        delete p;
    }
}

/**
 * Marker distributor, main goal is too evenly distribute the markers on the surface.
 * @bug marker dispersion seems uneven at times, investigate
 * @note This method is purely for ibm markers.
 * @attention surface is blind of its own inner structure
 */
void markerPoints::distribute_markers() {
    // start at the first surface and look
    if(sg != nullptr) {
        // find the length of the ibm surfaces
        double ibm_surface = sg->calculate_total_surface_length(IBM);
        double total = sg->calculate_total_surface_length();
        // the markers should be equally space over the total surface
        // we dont want any rest we just want equal spacing
        // determine the number of markers in a surface
        double markers_fit_in= std::floor(ibm_surface/marker_distance);
        // check weather or not we have to equalize so that everything is equally spaced
        // watch out for low marker numbers
        double add_on = std::fmod(ibm_surface,marker_distance);
        marker_distance += add_on/markers_fit_in;
        // we now step through the surface again
        // setup walk distance initial
        double walk_distance = marker_distance;
        for(auto s : sg->surfaces) {
            // only place markers on an ibm surface
            if(s->type != IBM) {
                continue;
            }
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

/**
 * Places marker on the periodics line.
 * @todo you buggy
 * @param line
 */
void markerPoints::distribute_markers_periodic(straight_t *line,boundaryType_t type, kernelType_t t) {
    // distributes markers for the use in periodic boundaries
    // used for placing of nodes and connectivity
    // pick the starter point and ceil it up
    // preliminary checks
    point_t startpoint = {std::ceil(line->point.x()),
                     std::ceil(line->point.y())};
    point_t helper = line->point + line->max_t*line->direction.normalized();
    point_t endpoint = {std::floor(helper.x()),
                        std::floor(helper.y())};
    // check the distance between start
    vector_t distance = endpoint - startpoint;
    if(distance.norm() < 1e-10) {
        std::cerr << "The straight line given is too short to work with." << std::endl;
        return;
    }
    // creat a new line with the right properties
    straight_t straight;
    // setup
    point_t current = startpoint;
    // manipulate the data based on kernel ibm
    if(type == IBM) {
        double ibm_range = kernel_id_to_lattice_search(t);
        // move the startpoint
        current = startpoint - ibm_range*line->direction;
        straight.point = current;
        straight.direction = line->direction;
        // move the end
        straight.max_t = distance.norm() + 2*ibm_range;
    }
    // now lets walk the walk
    vector_t walker = straight.direction.normalized();
    bool end = false;
    double dk = 0;
    while(!end) {
        if(point_on_straight(&straight,&current,&dk)) {
            // add a marker point
            auto add_me = new point_t;
            *add_me = current;
            marker_points.push_back(add_me);
        }
        else {
            end = true;
        }
        current += walker;
    }

}

/**
 * Gives out the actual marker distance used.
 * @return
 */
double markerPoints::return_marker_distance() {
    return marker_distance;
}

/**
 * Writes out the marker points into a file or not.
 * @param write_file
 */
void markerPoints::write_out_markers(bool write_file) {
    std::ofstream out;
    out.open("markers");
    if(out.is_open()) {
        for(auto m : marker_points) {
            out << m->x() << " " << m->y() << std::endl;
        }
    }
}

/**
 * Constructor.
 * @param in
 * @param out
 * @param ori
 */
periodicBundles::periodicBundles(markerPoints *in, markerPoints *out, bool ori)
    : inlet(in), outlet(out), orientation(ori)
{
    // check sizes
    if(inlet->marker_points.size() != outlet->marker_points.size()) {
        std::cerr << "Sizes of the inlet and outlet dont agree!" << std::endl;
    }
}

/**
 * Just makes and association between the position of the inlet/outlets position in the saving arrays.
 */
void periodicBundles::connect() {
    // loop over inlet and outlet
    size_t size = inlet->marker_points.size();
    for(size_t i = 0; i < size; ++i) {
        // when orientation is true the order is flipped
        size_t second = i;
        if(orientation) {
            second = size - 1 -i;
        }
        associations.emplace(i,second);
    }
}
