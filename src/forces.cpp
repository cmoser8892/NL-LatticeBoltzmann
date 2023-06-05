#include "forces.h"

/**
 * @fn inline double calculate_truncation_force(vector_t* c, vector_t* u, vector_t* force)
 * @brief calculates the ugly part of  th truncation of the forcoe term term (viggen 236)
 * @param c veloicty set
 * @param u velocitay
 * @param force
 * @return
 */
double calculate_truncation_force(vector_t* c, vector_t* u, vector_t* force) {
   double return_value = 0;
   double cs_2 = 1.0/3;
   for(int alpha = 0; alpha < c->size(); ++alpha) {
       for(int beta = 0; beta < u->size(); ++beta) {
           return_value += ((c->operator[](alpha)/cs_2) + (((c->operator[](alpha)*c->operator[](beta) - cs_2 * conical_delta(alpha,beta)) *u->operator[](beta)) / (cs_2*cs_2)))*force->operator[](alpha);
       }
   }
   // std::cout << std::endl;
   return return_value;
}

/// helper and sub-classes
/**
* @fn double gladrowForce::return_position_x(point_t *self_position)
* @brief
* @param self_position
* @return
*/
double gladrowForce::return_position_x(point_t *self_position) {
   // we determine the force magnitude
   vector_t distance_vector = (*self_position) - middle;
   double force_magnitude = distance_vector.norm()/max_distance * force;
   // we determine the x component
   distance_vector /= distance_vector.norm();
   // change directions and multiply with the force magnitude
   return distance_vector.y() * force_magnitude;
}

/**
* @fn double gladrowForce::return_position_y(point_t *self_position)
* @brief returns the force at a specific point
* @param self_position
* @return
*/
double gladrowForce::return_position_y(point_t *self_position) {
   // we determine the force magnitude
   vector_t distance_vector = (*self_position) - middle;
   double force_magnitude = distance_vector.norm()/max_distance * force;
   // we determine the x component
   distance_vector /= distance_vector.norm();
   // change directions and multiply with the force magnitude
   return - distance_vector.x() * force_magnitude;
}

/**
* @fn gladrowForce::gladrowForce(double f, point_t canvas_size)
* @brief constructor sets the size and middle of the rotating frame, also sets the basic force
* @param f
* @param canvas_size
*/
gladrowForce::gladrowForce(double f, point_t canvas_size) {
   // force denotes the force at the extremes of the canvas
   max_distance = (canvas_size - canvas_size/2).norm();
   middle = canvas_size/2;
   force = f;
}

/**
* @fn double gladrowForce::return_current_next_x(point_t *self_position, int channel)
* @brief function to return the next x foce component
* @param self_position
* @param channel
* @return
*/
double gladrowForce::return_current_next_x(point_t *self_position, int channel) {
   // setup the cidt position
   array_t intermediate = *self_position;
   point_t next = intermediate + velocity_set.col(channel);
   // calculate
   return return_position_x(self_position) + return_position_x(&next);
}

/**
* @fn double gladrowForce::return_current_next_y(point_t *self_position, int channel )
* @brief calculates the y forces
* @param self_position
* @param channel
* @return
*/
double gladrowForce::return_current_next_y(point_t *self_position, int channel ) {
   // setup the cidt position
   array_t intermediate = *self_position;
   point_t next = intermediate + velocity_set.col(channel);
   // calculate
   return return_position_y(self_position) + return_position_y(&next);
}

/// class goa force
/**
 * @fn inline double goaForce::truncation_force()
 * @brief non optimized version of the channel force calculation
 * @return
 */
inline double goaForce::truncation_force() {
   double return_value = 0;
   double cs_2 = 1.0/3;
   for(int alpha = 0; alpha < velocity_channel_set.size(); ++alpha) {
       for(int beta = 0; beta < velocity.size(); ++beta) {
           return_value += ((velocity_channel_set(alpha)/cs_2) + (((velocity_channel_set(alpha)*velocity_channel_set(beta) - cs_2 * (alpha == beta)) *velocity(beta)) / (cs_2*cs_2)))*force_alpha(alpha);
       }
   }
   // std::cout << std::endl;
   return return_value;
}

/**
 * @fn inline array_t goaForce::truncation_force_array()
 * @brief inlined the complicated part of the "channel" dependent force for 2dq9
 * @return
 */
inline array_t goaForce::truncation_force_array() {
   //
   array_t return_array;
   return_array.resize(9);
   /// channel 0
   return_array[0] = -3*velocity[0]*force_alpha[0] - 3*velocity[1]*force_alpha[1];
   /// channel 1 - 4
   return_array[1] =  6*force_alpha[0] + 6*velocity[0]*force_alpha[0] - 3*velocity[1]*force_alpha[1];
   return_array[2] = -3*velocity[0]*force_alpha[0] + 6*force_alpha[1]  + 6*velocity[1]*force_alpha[1];
   return_array[3] = -6*force_alpha[0] + 6*velocity[0]*force_alpha[0]  - 3*velocity[1]*force_alpha[1];
   return_array[4] = -3*velocity[0]*force_alpha[0] - 6*force_alpha[1] + 6*velocity[1]*force_alpha[1];
   /// channel 5 - 6
   return_array[5] =   6* force_alpha[0] + 6*velocity[0]*force_alpha[0] + 9*velocity[1]*force_alpha[0]
                     + 6* force_alpha[1] + 9*velocity[0]*force_alpha[1] + 6*velocity[1]*force_alpha[1];
   return_array[6] = - 6* force_alpha[0] + 6*velocity[0]*force_alpha[0] - 9*velocity[1]*force_alpha[0]
                     + 6* force_alpha[1] - 9*velocity[0]*force_alpha[1] + 6*velocity[1]*force_alpha[1];
   /// channel 7 - 8
   return_array[7] = - 6* force_alpha[0] + 6*velocity[0]*force_alpha[0] + 9*velocity[1]*force_alpha[0]
                     - 6* force_alpha[1] + 9*velocity[0]*force_alpha[1] + 6*velocity[1]*force_alpha[1];
   return_array[8] =   6* force_alpha[0] + 6*velocity[0]*force_alpha[0] - 9*velocity[1]*force_alpha[0]
                     - 6* force_alpha[1] - 9*velocity[0]*force_alpha[1] + 6*velocity[1]*force_alpha[1];;
   // return the array
   return return_array;
}

/**
 * @fn void goaForce::calculate_F_circle(point_t* p)
 * @brief circular force field
 * @param p
 */
void goaForce::calculate_F_circle(point_t* p, double max_force_magnitude,double ux, double uy) {
   velocity.x() = ux;
   velocity.y() = uy;
   double max_distance = (size - size/2).norm();;
   // we determine the force magnitude
   vector_t distance_vector = (*p) - middle;
   double force_magnitude = distance_vector.norm()/max_distance * max_force_magnitude;
   // we determine the x component
   distance_vector /= distance_vector.norm();
   // change directions and multiply with the force magnitude
   force_alpha = {  distance_vector.y() * force_magnitude,
                  - distance_vector.x() * force_magnitude};
}

/**
* @fn void goaForce::calculate_F_rotation()
* @brief calculates the physical force in the domain
*/
void goaForce::calculate_F_rotation(double ux, double uy, point_t* p) {
   // todo what to do with mass?!
   // setup the internal variables
   velocity.x() = ux;
   velocity.y() = uy;
   radius = calculate_distance(p,&origin);
   vector_t origin_to_point = *p - origin;
   angle = calculate_angle(&reference,&origin_to_point);
   // a_z =  w2r
   double force_z_alpha_value = omega * omega* radius;
   force_alpha.x() = 0;
   force_alpha.y() = force_z_alpha_value;
   // we rotate the force so its always out looking
   Eigen::Rotation2D<double> rot;
   rot.angle() = angle; // parts of pi use EIGEN_PI
   force_alpha = rot * force_alpha;
   // a_c = -2 (w x v)
    vector3d_t omega_vec = {0,0,omega};
    vector3d_t velocity_vec = {velocity.x(),velocity.y(),0};
    vector3d_t result = -2 *(velocity_vec.cross(omega_vec));
    // crash in case we get a z value
    assert(result.z() == 0);
    force_alpha.x() += result.x();
    force_alpha.y() += result.y();
}

/**
* @fn calculates the channel dependent force
* @brief
*/
void goaForce::calculate_F_i() {
   force_channels = truncation_force_array();
}

/**
* @fn goaForce::goaForce(point_t o, point_t c, double o1, double o2)
* @brief construtor does constructor things
* @param o
* @param c
* @param o1
* @param o2
*/
goaForce::goaForce(point_t o, point_t c, double o1) {
   origin = o;
   size = c;
   middle = c/2;
   omega = o1;
   force_channels.resize(CHANNELS);
   force_channels.setZero();
   reference ={1,0};
}

/**
* @fn vector_t goaForce::return_force_alpha()
* @brief returns the alpha force component
* @return
*/
vector_t goaForce::return_force_alpha() {
   // useful for testing
   return force_alpha;
}

/**
 * @fn void goaForce::set_force_alpha(vector_t f)
 * @brief sets the force value in the force class (for testing)
 * @param f
 */
void goaForce::set_force_alpha(vector_t f) {
   force_alpha = f;
}

/**
 * @fn void goaForce::set_velocity(vector_t v)
 * @brief sets tje velocity value in the force class (for testing)
 * @param v
 */
void goaForce::set_velocity(vector_t v) {
   velocity = v;
}

/**
 * @fn void goaForce::set_omega(double o)
 * @brief sets the omega rotation in plane
 * @param o
 */
void goaForce::set_omega(double o) {
   omega = o;
}