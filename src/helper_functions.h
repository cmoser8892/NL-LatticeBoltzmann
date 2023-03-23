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
// reduce 32 bit numbers to 2 bit numbers
uint32_t reduce_32_2(uint32_t);
/// helper classes/sub-classes
// watchdog for rho
class rhoWatchdog {
    /// just prints out std error msgs
  private:
    flowfield_t rho;
    double sensitivity = 0.1; // displacement to the previous value
  public:
    rhoWatchdog(double s,point_t size);
    bool check(node* n,int step);
};

class pointKeyHash {
  public:
    std::unordered_multimap<handle_t,handle_t> keys;
    void fill_key(handle_t positions_handle, point_t pos);
    handle_t key_translation(point_t pos);
    handle_t key_translation(coordinate_t cord);
};

#endif // NL_LATTICEBOLTZMANN_HELPER_FUNCTIONS_H
