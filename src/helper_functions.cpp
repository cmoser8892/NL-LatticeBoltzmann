
#include <fstream>
#include <iostream>
#include "helper_functions.h"

// compares two arrays point by point
bool node_position_comparison(node* n, array_t* position) {
    array_t* node_position = &n->position;
    // check same size
    if(node_position->size() != position->size()) {
        return false;
    }
    // check same values
    for(int i = 0; i < node_position->size(); ++i) {
        if(node_position->operator()(i) == position->operator()(i)) {}
        else {
            return false;
        }
    }
    return true;
}

// writes the data to
void write_flowfield_data(flowfield_t * field, std::string filename, bool write_to_file) {
    std::ofstream out;
    out.open(filename);
    // preamble size infos
    // print out the rows first then to file
    for(int i = 0; i < field->rows(); ++i) {
        if(write_to_file) {
            out << field->row(i) << std::endl;
        }
        else {
            std::cout << field->row(i) << std::endl;
        }
    }
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

bool same_index(node* n) {
    // array wise comparison
    array_t pos = n->position;
    for(auto n: pos) {
        for(auto comp : pos) {
            // == didnt know is a buildin c++ comparator
            if(n == comp) {
            }
            else {
                return false;
            }
        }
    }
    return true;
}