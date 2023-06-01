#ifndef NL_LATTICEBOLTZMANN_FORCES_H
#define NL_LATTICEBOLTZMANN_FORCES_H

#include "types.h"
#include "node.h"
#include "helper_functions.h"

// create circular force
std::vector<vector_t> circular_force_generation(int total_steps, int switch_time, double magnitude);
double calculate_truncation_force(array_t c, array_t u, vector_t force);

// classes
class circularForce {
    /*
     * also not used right now
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

class gladrowForce {
  private:
    double force = 0;
    point_t middle;
    double max_distance;
    double return_position_x(point_t* self_position);
    double return_position_y(point_t* self_position);
  public:
    // first implementation assumes that we rotate around the center of the canvas
    gladrowForce(double force, point_t canvas_size);
    double return_current_next_x(point_t* self_position, int channel);
    double return_current_next_y(point_t* self_position, int channel );
};

class goaForce {
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
    void calculate_F_circle(point_t* p);
    void calculate_F_alpha();
    void calculate_F_i();
  public:
    goaForce(point_t origin, point_t canvas_size, double omega_1, double omega_2);
    void precalculate(double ux, double uy,point_t* position);
    double return_force(int channel_i);
    vector_t return_force_alpha();
};


#endif // NL_LATTICEBOLTZMANN_FORCES_H
