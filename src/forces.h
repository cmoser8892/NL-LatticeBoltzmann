#ifndef NL_LATTICEBOLTZMANN_FORCES_H
#define NL_LATTICEBOLTZMANN_FORCES_H

#include "types.h"
#include "node.h"
#include "helper_functions.h"

// create circular force
double calculate_truncation_force(vector_t* c, vector_t* u, vector_t* force);

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
    vector_t velocity;
    vector_t force_alpha;
    vector_t velocity_channel_set;
    vector_t reference;
    // vars
    double radius = 0;
    double angle = 0;
    double omega = 0;
    // functions
    inline double truncation_force();
    inline array_t truncation_force_array();
  public:
    array_t force_channels;
    goaForce(point_t origin, point_t canvas_size, double omega);
    void calculate_F_circle(point_t* p, double max_force_magnitude, double ux, double uy);
    void calculate_F_rotation(double ux, double uy, point_t * p);
    void calculate_F_i();
    vector_t return_force_alpha();
    void set_force_alpha(vector_t f);
    void set_velocity(vector_t v);
    void set_omega(double o);
};


#endif // NL_LATTICEBOLTZMANN_FORCES_H
