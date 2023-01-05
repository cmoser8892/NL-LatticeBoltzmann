//
// Created by workstation on 16.11.22.
//

#ifndef MY_LATTICE_BOLTZMANN_H
#define MY_LATTICE_BOLTZMANN_H

// includes
#include "boundary_point_generator.h"
#include "node.h"
#include "nodeGenerator.h"
#include "types.h"

// for completeness
extern matrix_t velocity_set;
// functions
typedef array_t (*equilibrium_function) (node* node);
typedef void (*streaming_function) (node* node); // could also do boundary conditions there?!
typedef void (*collision_function) (node* node);
typedef void (*calculate_macro_values) (node* node);

class simulation {
  private:
    double relaxation;
    double u_wall;
    boundaryPointConstructor * boundary_points = nullptr;
    nodeGenerator* node_generator = nullptr;
    void stream_links(node* n);
    int links_correct_channel(node* n, int link_channel);
    int switch_link_dimensions(int link_channel);
    double extra_moving_bb_equation(); // todo me now
  public:
    std::vector<node*> nodes;
    simulation(boundaryPointConstructor* c);
    simulation(boundaryPointConstructor* c,nodeGenerator* g);
    void streaming_step_1();
    void bounce_back();
    void streaming_step_2();
    void init();
    void run();
    void get_data(bool write_to_file);
    std::vector<node*> access();
};

#endif // MY_LATTICE_BOLTZMANN_H
