//
// Created by workstation on 16.11.22.
//

#ifndef MY_LATTICE_BOLTZMANN_H
#define MY_LATTICE_BOLTZMANN_H

// includes
#include "types.h"
#include "node.h"
#include "boundary_point_generator.h"
#include "node_generator.h"

// for completeness
extern matrix_t velocity_set;
// functions
typedef array_t (*equilibrium_function) (node* node);
typedef void (*streaming_function) (node* node); // could also do boundary conditions there?!
typedef void (*collision_function) (node* node);
typedef void (*calculate_macro_values) (node* node);

class simulation {
  private:
    std::vector<node*> nodes;
    int dimensions, channels;
    int size_x,size_y; // prob best to put them into arrays
    int limit_x, limit_y;
  public:
    simulation() = default;
    void init(int six, int siy);
    void run();
    void get_data();
    std::vector<node*> access();
};

#endif // MY_LATTICE_BOLTZMANN_H
