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

typedef struct simulation_parameters {
    double relaxation;
    double u_wall;
}simulation_parameters_t;

class simulation {
  private:
    double relaxation;
    double u_wall;
    boundaryPointConstructor * boundary_points = nullptr;
    nodeGenerator* node_generator = nullptr;
  public:
    std::vector<node*> nodes;
    explicit simulation(boundaryPointConstructor* c);
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
