//
// Created by simulation on 21.11.22.
//

#ifndef NL_LATTICEBOLTZMANN_FUNCTIONS_H
#define NL_LATTICEBOLTZMANN_FUNCTIONS_H

#include "types.h"
#include "lattice_boltzmann.h"

array_t equilibrium(node* node);
void streaming(node* node);
void collision(node* node);
void macro(node * node);
#endif // NL_LATTICEBOLTZMANN_FUNCTIONS_H
