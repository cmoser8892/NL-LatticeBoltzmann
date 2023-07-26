#ifndef NL_LATTICEBOLTZMANN_FUNCTIONS_H
#define NL_LATTICEBOLTZMANN_FUNCTIONS_H

#include "types.h"
#include "node.h"
// functions is specific for nodes
// functional methods
array_t equilibrium(node* node);
array_t equilibrium_general(node* node);
double calculate_later_equilibrium(double cx, double cy, double ux, double uy);
array_t equilibrium_2d(double ux, double uy, double rho);
void collision(node* node,double relaxation);
void fused_collision(node* node, double relaxation);
void macro(node * node);
void fused_macro(node* node);
// todo still being worked on
void pressure_periodic_in(node* node, double rho);
void pressure_periodic_out(node* node, double rho);
double bb_switch_channel(int from_channel, double uw);
int switch_link_dimensions(int link_channel);
std::tuple<double,double,double> calculate_force_macro_values(array_t a, array_t f);
// writeing methods
void write_rho(node* node, flowfield_t* rho);
void write_ux(node* node, flowfield_t* uy); // writes velocity based on the place in the field
void write_uy(node* node,flowfield_t* ux);
double calculate_rho(oNode* n, int offset);
// unified calculate method
std::tuple<double,double,double> calculate_macro_population(array_t* p);
std::tuple<double, double> test_calculate_rho_both(array_t* p);
// kernel functions
double kernel_A(double range);
double kernel_B(double range);
double kernel_C(double range);
double kernel_A_2d(vector_t* p);
double kernel_B_2d(vector_t* p);
double kernel_C_2d(vector_t *p);
double kernel_id_to_lattice_search(kernelType_t t);
bool point_on_straight(straight_t * s, point_t *p, double* overshot);
vector_t compute_lagrangian_force();
// debug method
void debug_node(node* node, bool printing);
#endif // NL_LATTICEBOLTZMANN_FUNCTIONS_H