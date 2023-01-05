//
// Created by simulation on 21.11.22.
//

#ifndef NL_LATTICEBOLTZMANN_FUNCTIONS_H
#define NL_LATTICEBOLTZMANN_FUNCTIONS_H

#include "types.h"
#include "node.h"

extern double relaxation;
array_t equilibrium(node* node);
void collision(node* node,double relaxation);
void macro(node * node);
double bb_switch_channel(int link_channel, double uw);
// writeing methods
void write_rho(node* node, flowfield_t* rho);
void write_ux(node* node, flowfield_t* uy); // writes velocity based on the place in the field
void write_uy(node* node,flowfield_t* ux);
// debug method
void debug_node(node* node, bool printing);
#endif // NL_LATTICEBOLTZMANN_FUNCTIONS_H