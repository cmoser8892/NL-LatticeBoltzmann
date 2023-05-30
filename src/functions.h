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
void pressure_periodic_in(node* node, double rho);
void pressure_periodic_out(node* node, double rho);
double bb_switch_channel(int from_channel, double uw);
int switch_link_dimensions(int link_channel);
// writeing methods
void write_rho(node* node, flowfield_t* rho);
void write_ux(node* node, flowfield_t* uy); // writes velocity based on the place in the field
void write_uy(node* node,flowfield_t* ux);
double calculate_ux(oNode* n,int offset);
double calculate_uy(oNode* n, int offset);
double calculate_rho(oNode* n, int offset);
// debug method
void debug_node(node* node, bool printing);
#endif // NL_LATTICEBOLTZMANN_FUNCTIONS_H