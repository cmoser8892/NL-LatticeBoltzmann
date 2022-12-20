
#include <fstream>
#include <iostream>
#include "helper_functions.h"

// compares two arrays point by point
bool compare_arrays(array_t a1, array_t a2) {
    if(a1.size() != a2.size()){
        return false;
    }
    else {
        // comparing doubles is cringe
        for(int i = 0; i < a1.size(); ++i) {
            if(std::abs(a1(i)- a2(i)) < 1.0/DOUBLE_COMPARISON_PRECISION) {}
            else {
                return false;
            }
        }
    }
    return true;
}

// writes the data to
void write_flowfield_data(flowfield_t * field, std::string filename) {

    std::ofstream out;
    out.open(filename);
    // preamble size infos
    // print out the rows first then to file
    for(int i = 0; i < field->rows(); ++i)
        std::cout << field->row(i) << std::endl;
    std::cout << std::endl;
    out.close();
}

bool check_inside_limits_upper_lower(point_t* p, point_t* limit_lower, point_t* limit_upper) {
    bool return_value = true;
    if(p->x() < limit_lower->x() || p->x() > limit_upper->x()) {
        return_value = false;
    }
    if(p->y() < limit_lower->y() || p->y() > limit_upper->y()) {
        return_value = false;
    }
    return return_value;
}

bool compare_two_points(point_t* p1, point_t* p2 ) {
    // todo double comparison :)(
    bool return_value = true;
    if( p1->x() != p2->x()) {
        return_value = false;
    }
    if( p1->y() != p2->y()) {
        return_value = false;
    }
    return return_value;
}