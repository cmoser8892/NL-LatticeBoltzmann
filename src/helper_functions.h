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
double calculate_distance(point_t* p1, point_t* p2);
double calculate_truncation_force(array_t c, array_t u, vector_t force);
int conical_delta(int a, int b);
// base path
std::filesystem::path get_base_path();
// create circular force
std::vector<vector_t> circular_force_generation(int total_steps, int switch_time, double magnitude);
/// helper classes/sub-classes
// circular force as a class instead of a big vector
class circularForce {
    /*
     * This force hides the housekeeping necessary for the calculation of a force
     * going around and switching directions ever so often, it is supposed to be just
     * an adhoc test to get a familiar with the collision term with a force
     */
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
    double return_current_next_x();
    double return_current_next_y();
    void increment();
};

class circleForce {
  private:
    double force = 0;
    point_t middle;
    double max_distance;
    double return_position_x(point_t* self_position);
    double return_position_y(point_t* self_position);
  public:
    // first implementation assumes that we rotate around the center of the canvas
    circleForce(double force, point_t canvas_size);
    double return_current_next_x(point_t* self_position, int channel);
    double return_current_next_y(point_t* self_position, int channel );
};

class rotatingForce {
  private:
    point_t origin;
    point_t size;
    point_t middle;
    vector_t omega;
    vector_t velocity;
    vector_t force_alpha;

    double radius = 0;
    double angle = 0;
    // 9 long for 2dq9
    array_t force_channels;
    void calculate_F_alpha();
    void calculate_F_i();
  public:
    rotatingForce(point_t origin, point_t canvas_size, double omega_1, double omega_2);
    void precalculate(double ux, double uy,point_t* position);
    double return_force(int channel_i);
    vector_t return_force_alpha();
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
