//
// Created by christoph on 23.11.22.
//

#ifndef NL_LATTICEBOLTZMANN_HELPER_FUNCTIONS_H
#define NL_LATTICEBOLTZMANN_HELPER_FUNCTIONS_H

#include "types.h"
#include "node.h"
#include <iostream>
#include <sstream>
#include <vector>

bool node_position_comparison(node* n, array_t* pos);
void write_flowfield_data(flowfield_t* field,std::string filename,bool write_to_file);
bool check_inside_limits_upper_lower(point_t* p, point_t* limit_lower, point_t* limit_upper);
bool compare_two_points(point_t* p1, point_t* p2 );
bool same_index(node* n);
std::vector<std::string> split_string (std::string s, std::string delimiter);
// bit interleaving
uint64_t bit_interleaving(uint32_t x, uint32_t y);
uint64_t bit_interleaving_2d(uint32_t, uint32_t);
uint64_t bit_interleaving_3d(uint32_t, uint32_t, uint32_t);
// bit extraleaving
uint32_t bit_extraleaving_2d_x(uint64_t);
uint32_t bit_extraleaving_2d_y(uint64_t);
uint32_t bit_extraleaving_3d_x(uint64_t);
uint32_t bit_extraleaving_3d_y(uint64_t);
uint32_t bit_extraleaving_3d_z(uint64_t);
// handle compare for sorting
bool compare_handles(nodePoint_t * a, nodePoint_t *b);
#endif // NL_LATTICEBOLTZMANN_HELPER_FUNCTIONS_H
