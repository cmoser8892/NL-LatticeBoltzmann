#include "forces.h"

/**
 * @fn double calculate_truncation_force(array_t c, array_t u, vector_t force)
 * @brief calculates the ugly part of  th truncation of the forcoe term term (viggen 236)
 * @param c veloicty set
 * @param u velocitay
 * @param force
 * @return
 */
double calculate_truncation_force(array_t c, array_t u, vector_t force) {
   double return_value = 0;
   double cs_2 = 1.0/3;
   for(int alpha = 0; alpha < c.size(); ++alpha) {
       for(int beta = 0; beta < u.size(); ++beta) {
           // std::cout << c(alpha) << std::endl;
           // std::cout << c(beta) << std::endl;
           return_value += ((c(alpha)/cs_2) +
                            (((c(alpha)*c(beta) - cs_2 * conical_delta(alpha,beta)) *u(beta)) / (cs_2*cs_2)))
                           *force(alpha);
           // std::cout << return_value << std::endl;
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

/// class rotating force
void goaForce::calculate_F_circle(point_t* p) {
   double max_distance = (size - size/2).norm();;
   double force = 0.0035;
   // we determine the force magnitude
   vector_t distance_vector = (*p) - middle;
   double force_magnitude = distance_vector.norm()/max_distance * force;
   // we determine the x component
   distance_vector /= distance_vector.norm();
   // change directions and multiply with the force magnitude
   force_alpha = {  distance_vector.y() * force_magnitude,
                  - distance_vector.x() * force_magnitude};
}

/**
* @fn void goaForce::calculate_F_alpha()
* @brief calculates the physical force in the domain
*/
void goaForce::calculate_F_alpha() {
   force_alpha.setZero();
   // F_c = -2 (w x v)
   double force_c_alpha_value = -2 * omega.norm() * velocity.norm();
   double force_z_alpha_value = omega.norm() * omega.norm() * radius;
   // F_z =  w2r
   force_alpha.x() = force_c_alpha_value;
   force_alpha.y() = force_z_alpha_value;
   // F_c = -2 (w x v)
   Eigen::Rotation2D<double> rot;
   rot.angle() = angle; // parts of pi use EIGEN_PI
   force_alpha = rot * force_alpha;
}

/**
* @fn calculates the channel dependent force
* @brief
*/
void goaForce::calculate_F_i() {
   for(int i = 0; i < CHANNELS; ++i) {
       force_channels[i] = weights(i) * calculate_truncation_force(velocity_set.col(i),velocity,force_alpha);
   }
}

/**
* @fn goaForce::goaForce(point_t o, point_t c, double o1, double o2)
* @brief construtor does constructor things
* @param o
* @param c
* @param o1
* @param o2
*/
goaForce::goaForce(point_t o, point_t c, double o1, double o2) {
   // todo might be a good idea to reduce omega doesnt make sense as a vector
   origin = o;
   size = c;
   middle = c/2;
   omega.x() = o1;
   omega.y() = o2;
   force_channels.resize(CHANNELS);
   force_channels.setZero();
}

/**
* @fn void goaForce::precalculate(double ux, double uy, point_t *position)
* @brief precalculates everything in used in the fused macro streaming step
* @param ux
* @param uy
* @param position
*/
void goaForce::precalculate(double ux, double uy, point_t *position) {
   // set velocity
   velocity.x() = ux;
   velocity.y() = uy;
   // calculate distance
   // radius = calculate_distance(position,&origin);
   // angle = 2* asin(calculate_distance(&middle,position)/(2*radius));
   calculate_F_circle(position);
   // calculate_F_alpha();
   calculate_F_i();
}

/**
* @fn double goaForce::return_force(int channel_i)
* @brief calculates the forces dependent on the channels
* @param channel_i
* @return
*/
double goaForce::return_force(int channel_i) {
   // truncation of the forcing term (viggen 236)
   return force_channels[channel_i];
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