//
// Created by workstation on 16.11.22.
//

#ifndef MY_LATTICE_BOLTZMANN_H
#define MY_LATTICE_BOLTZMANN_H

// includes
#include "types.h"
// for completeness
class node;
extern array_t velocity_set;
// functions
typedef array_t (*equilibrium_function) (node* node);
typedef void (*streaming_function) (node* node); // could also do boundary conditions there?!
typedef void (*collision_function) (node* node);
typedef void (*calculate_macro_values) (node* node);

class node {
  public:
    node_identifier_t node_type;
    array_t data; // doesnt have a set size; maybe could also do with a vector
    array_t copy;
    std::vector<node*> neighbors;
    // macro values are local
    double rho;
    array_t u;
    array_t position;
    // pointer to functions no idea if they can work on the data and if i should
    // implemement that here?!
    equilibrium_function equilibrium_func;
    streaming_function streaming_func;
    collision_function collision_func;
    calculate_macro_values macro_func;
    // methods
    node(int dimensions, int channels, array_t positions);
};

class simulation {
  private:
    std::vector<node*> nodes;
    int dimensions, channels;
  public:
    simulation() = default;
    void init(void);
    void run(void);
};

#endif // MY_LATTICE_BOLTZMANN_H
