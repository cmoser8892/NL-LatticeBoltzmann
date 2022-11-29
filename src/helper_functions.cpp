
#include <fstream>
#include <iostream>
#include "helper_functions.h"

// compares two arrays point by point
bool compare_arrays(array_t a1, array_t a2) {
    if(a1.size() != a2.size()){
        return false;
    }
    else {
        for(int i = 0; i < a1.size(); ++i) {
            if(a1(i) != a2(i)) {
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
        out << field->row(i) << std::endl;
    out.close();
}