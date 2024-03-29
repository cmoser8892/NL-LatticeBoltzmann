#ifndef NL_LATTICEBOLTZMANN_HELPER_FUNCTIONS_H
#define NL_LATTICEBOLTZMANN_HELPER_FUNCTIONS_H

#include "types.h"
#include "node.h"
#include <iostream>
#include <sstream>
#include <vector>
#include <list>
#include <filesystem>

//vector nonsense
bool node_position_comparison(node* n, array_t* pos);
void write_flowfield_data(flowfield_t* field,std::string filename,bool write_to_file);
bool check_inside_limits_upper_lower(point_t* p, point_t* limit_lower, point_t* limit_upper);
bool compare_two_points(point_t* p1, point_t* p2 );
std::vector<std::string> split_string (std::string s, std::string delimiter);
bool check_plus_minus_90(vector_t * to_be_checked, vector_t* reference);
vector_t vector_to_cardinal(vector_t reference);
int index_of_velocity_set(vector_t set);
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
// reduce 32 bit numbers to 2 bit numbers
uint32_t reduce_32_2(uint32_t);
double calculate_distance(point_t* p1, point_t* p2);
int conical_delta(int a, int b);
double calculate_angle(vector_t* v1, vector_t* v2);
coordinate_t floor_position(point_t position);
coordinate_t add_coordinates(coordinate_t &a, coordinate &b);
uint8_t point_on_boarder(point_t* p, vector_t* canvas_size);
// base path
std::filesystem::path get_base_path();
// solve
void solve_truncation_force_symbols(int channel);
#endif // NL_LATTICEBOLTZMANN_HELPER_FUNCTIONS_H
