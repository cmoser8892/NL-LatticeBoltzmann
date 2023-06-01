#include "forces.h"

/**
* @fn std::vector<vector_t> circular_force_generation(int total_steps, int switch_time, double magnitude)
* @brief
* @param total_steps
* @param switch_time
* @param magnitude
* @return
*/
std::vector<vector_t> circular_force_generation(int total_steps, int switch_time, double magnitude) {
   //
   std::vector<vector_t> return_vector;
   matrix_t helper = {{1,0,-1,0},
                      {0,1,0,-1}};
   int selector = 0;
   for(int i = 0; i < total_steps; ++i) {
       // switch selector up and down
       if((total_steps % switch_time) == 0) {
           selector = (++selector) % 3;
       }
       // we select on of the vectors and multiply with the magnitude
       vector_t current = helper.col(selector) * magnitude;
       // push back
       return_vector.push_back(current);
   }
   return return_vector;
}

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
* @fn circularForce::circularForce(long switchtime, double mag)
* @brief constructor of the circular force calculation
* @param switchtime
* @param mag
*/
circularForce::circularForce(long switchtime, double mag) {
   switch_time = switchtime;
   magnitude = mag;
   counter = 1;
   // setup x and y force
   x_force = {{1,0,-1,0}};
   y_force = {{0,1,0,-1}};
   x_force *= magnitude;
   y_force *= magnitude;
   // setup selectors
   selector_switchero();
}

/**
* @fn void circularForce::selector_switchero()
* @brief count up to three then go back to 0
*/
void circularForce::selector_switchero() {
   if((counter % switch_time) == 0) {
       current_selector =  (++current_selector)%3;
   }
   if(((counter +1) % switch_time) == 0) {
       next_selector =  (++next_selector)%3;
   }
}

/**
* @fn double circularForce::return_current_x()
* @brief returns th current x force component
* @return
*/
double circularForce::return_current_x() {
   return x_force(current_selector);
}

/**
* @fn double circularForce::return_current_y()
* @brief returns the current y force component
* @return
*/
double circularForce::return_current_y() {
   return y_force(current_selector);
}

/**
* @fn double circularForce::return_next_x()
* @brief returns the next x force component
* @return
*/
double circularForce::return_next_x() {
   return x_force(next_selector);
}

/**
* @fn double circularForce::return_next_y()
* @brief returns the next y force component
* @return
*/
double circularForce::return_next_y() {
   return y_force(next_selector);
}

/**
* @fn double circularForce::return_current_next_x()
* @brief returns both force components
* @return
*/
double circularForce::return_current_next_x() {
   return x_force(current_selector) + x_force(next_selector);
}

/**
* @fn double circularForce::return_current_next_y()
* @brief returns both force components
* @return
*/
double circularForce::return_current_next_y() {
   return y_force(current_selector) + y_force(next_selector);
}

/**
* @fn void circularForce::increment()
* @brief increments the counter and checks if we have to change the selectors
*/
void circularForce::increment() {
   // look if we need to change the selectors
   selector_switchero();
   ++counter;
}

/**
* @fn double circleForce::return_position_x(point_t *self_position)
* @brief
* @param self_position
* @return
*/
double circleForce::return_position_x(point_t *self_position) {
   // we determine the force magnitude
   vector_t distance_vector = (*self_position) - middle;
   double force_magnitude = distance_vector.norm()/max_distance * force;
   // we determine the x component
   distance_vector /= distance_vector.norm();
   // change directions and multiply with the force magnitude
   return distance_vector.y() * force_magnitude;
}

/**
* @fn double circleForce::return_position_y(point_t *self_position)
* @brief returns the force at a specific point
* @param self_position
* @return
*/
double circleForce::return_position_y(point_t *self_position) {
   // we determine the force magnitude
   vector_t distance_vector = (*self_position) - middle;
   double force_magnitude = distance_vector.norm()/max_distance * force;
   // we determine the x component
   distance_vector /= distance_vector.norm();
   // change directions and multiply with the force magnitude
   return - distance_vector.x() * force_magnitude;
}

/**
* @fn circleForce::circleForce(double f, point_t canvas_size)
* @brief constructor sets the size and middle of the rotating frame, also sets the basic force
* @param f
* @param canvas_size
*/
circleForce::circleForce(double f, point_t canvas_size) {
   // force denotes the force at the extremes of the canvas
   max_distance = (canvas_size - canvas_size/2).norm();
   middle = canvas_size/2;
   force = f;
}

/**
* @fn double circleForce::return_current_next_x(point_t *self_position, int channel)
* @brief function to return the next x foce component
* @param self_position
* @param channel
* @return
*/
double circleForce::return_current_next_x(point_t *self_position, int channel) {
   // setup the cidt position
   array_t intermediate = *self_position;
   point_t next = intermediate + velocity_set.col(channel);
   // calculate
   return return_position_x(self_position) + return_position_x(&next);
}

/**
* @fn double circleForce::return_current_next_y(point_t *self_position, int channel )
* @brief calculates the y forces
* @param self_position
* @param channel
* @return
*/
double circleForce::return_current_next_y(point_t *self_position, int channel ) {
   // setup the cidt position
   array_t intermediate = *self_position;
   point_t next = intermediate + velocity_set.col(channel);
   // calculate
   return return_position_y(self_position) + return_position_y(&next);
}

/// class rotating force
void rotatingForce::calculate_F_circle(point_t* p) {
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
* @fn void rotatingForce::calculate_F_alpha()
* @brief calculates the physical force in the domain
*/
void rotatingForce::calculate_F_alpha() {
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
void rotatingForce::calculate_F_i() {
   for(int i = 0; i < CHANNELS; ++i) {
       force_channels[i] = weights(i) * calculate_truncation_force(velocity_set.col(i),velocity,force_alpha);
   }
}

/**
* @fn rotatingForce::rotatingForce(point_t o, point_t c, double o1, double o2)
* @brief construtor does constructor things
* @param o
* @param c
* @param o1
* @param o2
*/
rotatingForce::rotatingForce(point_t o, point_t c, double o1, double o2) {
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
* @fn void rotatingForce::precalculate(double ux, double uy, point_t *position)
* @brief precalculates everything in used in the fused macro streaming step
* @param ux
* @param uy
* @param position
*/
void rotatingForce::precalculate(double ux, double uy, point_t *position) {
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
* @fn double rotatingForce::return_force(int channel_i)
* @brief calculates the forces dependent on the channels
* @param channel_i
* @return
*/
double rotatingForce::return_force(int channel_i) {
   // truncation of the forcing term (viggen 236)
   return force_channels[channel_i];
}

/**
* @fn vector_t rotatingForce::return_force_alpha()
* @brief returns the alpha force component
* @return
*/
vector_t rotatingForce::return_force_alpha() {
   // useful for testing
   return force_alpha;
}