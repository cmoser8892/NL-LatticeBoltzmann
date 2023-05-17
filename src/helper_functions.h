#ifndef NL_LATTICEBOLTZMANN_HELPER_FUNCTIONS_H
#define NL_LATTICEBOLTZMANN_HELPER_FUNCTIONS_H

#include "types.h"
#include "node.h"
#include <iostream>
#include <sstream>
#include <vector>
#include <list>
#include <filesystem>

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
// base path
std::filesystem::path get_base_path();
// create circular force
std::vector<vector_t> circular_force_generation(int total_steps, int switch_time, double magnitude);
/// helper classes/sub-classes
// circular force as a class instead of a big vector
class circularForce {
    // the whole purpuse of this class is to hide certain elements of the force
    // calculations
  private:
    long counter = 0;
    int current_selector =  0;
    int next_selector = 0;
    double magnitude = 0;
    long switch_time = 0;
    matrix_t x_force;
    matrix_t y_force;
    void selector_switchero();
  public:
    circularForce(long switchtime, double magnitude);
    double return_current_x();
    double return_current_y();
    double return_next_x();
    double return_next_y();
    void increment();
};

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

// works with bit interleaving floors points to a full point on the gird!
class pointKeyHash {
  private:
    std::unordered_multimap<handle_t,handle_t> keys;
  public:
    void fill_key(handle_t positions_handle, point_t pos);
    void clear();
    handle_t key_translation(point_t pos);
    handle_t key_translation(coordinate_t cord);
};

// basic stash that holds items for a number of iterations
// a window with a conveyor belt behind it, specifically for handles
class windowedHandles{
    // dont make this to big
  private:
    unsigned long target_size;
    std::list<handle_t> previous;
  public:
    // functions
    windowedHandles(unsigned long size);
    unsigned long size();
    void add(handle_t h);
    bool check(handle_t h);
};

#endif // NL_LATTICEBOLTZMANN_HELPER_FUNCTIONS_H
