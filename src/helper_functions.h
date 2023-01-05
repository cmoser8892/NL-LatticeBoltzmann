//
// Created by christoph on 23.11.22.
//

#ifndef NL_LATTICEBOLTZMANN_HELPER_FUNCTIONS_H
#define NL_LATTICEBOLTZMANN_HELPER_FUNCTIONS_H

#include "types.h"
#include "node.h"


bool node_position_comparison(node* n, array_t* pos);
void write_flowfield_data(flowfield_t* field,std::string filename,bool write_to_file);
bool check_inside_limits_upper_lower(point_t* p, point_t* limit_lower, point_t* limit_upper);
bool compare_two_points(point_t* p1, point_t* p2 );
bool same_index(node* n);
#endif // NL_LATTICEBOLTZMANN_HELPER_FUNCTIONS_H
