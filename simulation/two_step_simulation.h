#ifndef NL_LATTICEBOLTZMANN_TWO_STEP_SIMULATION_H
#define NL_LATTICEBOLTZMANN_TWO_STEP_SIMULATION_H

// includes
#include "boundary_point_generator.h"
#include "node.h"
#include "nodeGenerator.h"
#include "types.h"
#include "functions.h"
#include "forces.h"

// basic sim
class basicSimulation {
  private:
    simulation_parameters_t parameters;
    boundaryPointConstructor * boundary_points = nullptr;
    nodeGenerator* node_generator = nullptr;
  public:
    std::vector<node*> nodes;
    explicit basicSimulation(boundaryPointConstructor* c);
    ~basicSimulation();
    basicSimulation(boundaryPointConstructor* c,nodeGenerator* g);
    void set_simulation_parameters(simulation_parameters_t t);
    void streaming_step_1();
    void bounce_back();
    void streaming_step_2();
    void collisions();
    void fused_streaming(node* n);
    void fused_bounce_back(node* n);
    void init();
    void fused_init();
    void run();
    void fused_run();
    void get_data(bool write_to_file, point_t org);
    void delete_nodes();
};
#endif // NL_LATTICEBOLTZMANN_TWO_STEP_SIMULATION_H
