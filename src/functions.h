//
// Created by simulation on 21.11.22.
//

#ifndef NL_LATTICEBOLTZMANN_FUNCTIONS_H
#define NL_LATTICEBOLTZMANN_FUNCTIONS_H

#include "types.h"
#include "lattice_boltzmann.h"
extern double relaxation;
array_t equilibrium(node* node);
void streaming_step1(node* node);
void streaming_step2(node* node);
void collision(node* node);
void macro(node * node);
void moving_wall(node * node,int side_position,double uw);
// writeing methods
void write_ux(node* node, flowfield_t* uy); // writes velocity based on the place in the field
void write_uy(node* node,flowfield_t* ux);
// debug method
void debug_node_neighbors(node* node);
#endif // NL_LATTICEBOLTZMANN_FUNCTIONS_H