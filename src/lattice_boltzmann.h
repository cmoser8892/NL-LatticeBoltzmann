//
// Created by workstation on 16.11.22.
//

#ifndef MY_LATTICE_BOLTZMANN_H
#define MY_LATTICE_BOLTZMANN_H

// includes
#include "types.h"
// for completeness
class node;
typedef double (*equilibrium_function) (node* node);
typedef void (*streaming) (node* node); // could also do boundary conditions there?!
typedef void (*collision) (node* node);
typedef void (*calculate_macro_values) (node* node);

class node {
  private:
    node_identifier_t node_type;
    array_t data; // doesnt have a set size; maybe could also do with a vector
    std::vector<node*> neighbors;
    // macro values are local
    double rho;
    array_t u;
    array_t position;
    // pointer to functions no idea if they can work on the data and if i should
    // implemnte that here?!
    equilibrium_function equilibrium_function;
    streaming streaming_function;
    collision collision;
  public:
    node(int dimensions, int channels, array_t positions);
    void run();
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
