#ifndef NL_LATTICEBOLTZMANN_FORCES_H
#define NL_LATTICEBOLTZMANN_FORCES_H

#include "types.h"
#include "node.h"
#include "helper_functions.h"

// create circular force
double calculate_truncation_force(vector_t* c, vector_t* u, vector_t* force);
/**
 * Class the implements a circular force after Gladrow.
 * @ref DOI:10.1007/b72010
 */
class gladrowForce {
  private:
    double force = 0; /**< force value, gets set by the constructor */
    point_t middle; /**< middle point */
    double max_distance; /**< the force is maximal here (the set force value) */
    double return_position_x(point_t* self_position);
    double return_position_y(point_t* self_position);
  public:
    // first implementation assumes that we rotate around the center of the canvas
    gladrowForce(double force, point_t canvas_size);
    double return_current_next_x(point_t* self_position, int channel);
    double return_current_next_y(point_t* self_position, int channel );
};

/**
 * Class that implements parts (minus weights) of the force calculation
 * @ref PHYSICAL REVIEW E 66, 036304 (2002)
 */
class goaForce {
  private:
    point_t origin; /**<  origin of the rotation */
    point_t size; /**<  size of the simulation */
    point_t middle; /**<  middle point */
    vector_t velocity; /**<  velocity of the particle */
    vector_t force_alpha; /**<  2d force */
    vector_t velocity_channel_set; /**<  specific vector velocity */
    vector_t reference; /**<  reference to the angle calculation */
    // vars
    double radius = 0; /**< radius of our rotating frame */
    double angle = 0; /**< angle of from the frame */
    double omega = 0; /**< circular velocity */
    // functions
    inline void truncation_force(int channel);
    inline array_t truncation_force_array();
  public:
    array_t force_channels; /**<  truncation of the force (minus the weight) */
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
