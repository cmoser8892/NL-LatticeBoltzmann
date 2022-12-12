//
// Created by christoph on 23.11.22.
//

#ifndef NL_LATTICEBOLTZMANN_HELPER_FUNCTIONS_H
#define NL_LATTICEBOLTZMANN_HELPER_FUNCTIONS_H

#include "types.h"

// aka 1.0/100000 or whatever stands here
#define DOUBLE_COMPARISON_PRECISION 100000

// compare an array with a point
bool compare_arrays(array_t a1, array_t a2);
void write_flowfield_data(flowfield_t* field,std::string filename);

#endif // NL_LATTICEBOLTZMANN_HELPER_FUNCTIONS_H
