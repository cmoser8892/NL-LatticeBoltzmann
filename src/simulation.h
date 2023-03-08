#ifndef MY_LATTICE_BOLTZMANN_H
#define MY_LATTICE_BOLTZMANN_H

// includes
#include "boundary_point_generator.h"
#include "node.h"
#include "nodeGenerator.h"
#include "types.h"
#include "functions.h"

typedef struct simulation_parameters {
    double relaxation = 0.5;
    double u_wall = 0;
    bool bypass_lookup = true;
    double lookup_floor = -1.0;
    double lookup_ceiling = 1.0;
    uint32_t lookup_bits = 8;
}simulation_parameters_t;

class simulation {
  private:
    simulation_parameters_t parameters;
    boundaryPointConstructor * boundary_points = nullptr;
    nodeGenerator* node_generator = nullptr;
  public:
    std::vector<node*> nodes;
    explicit simulation(boundaryPointConstructor* c);
    ~simulation();
    simulation(boundaryPointConstructor* c,nodeGenerator* g);
    void set_simulation_parameters(simulation_parameters_t t);
    void streaming_step_1();
    void bounce_back();
    void streaming_step_2();
    void collisions();
    void init();
    void run();
    void fused_run();
    void get_data(bool write_to_file, point_t org);
    void delete_nodes();
};

#endif // MY_LATTICE_BOLTZMANN_H
