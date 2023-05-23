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
    void fused_streaming(node* n);
    void fused_bounce_back(node* n);
    void init();
    void fused_init();
    void run();
    void fused_run();
    void get_data(bool write_to_file, point_t org);
    void delete_nodes();
};

class oSimu {
  private:
    simulation_parameters_t parameters;
    boundaryPointConstructor * boundary_points = nullptr;
    nodeGenerator* node_generator = nullptr;
    circleForce* force;
  public:
    int offset_sim = 1;
    std::vector<oNode*> nodes;
    oSimu(boundaryPointConstructor* c,nodeGenerator* g);
    ~oSimu();
    void set_simulation_parameters(simulation_parameters_t t);
    void streaming(oNode* n);
    void bounce_back_moving(oNode* n);
    void one_step_macro_collision_forcing(oNode* n);
    void init();
    void run(int current_step);
    void current_run(int current_step);
    void get_data(bool write_to_file, point_t org);
    void get_data(bool write_to_file);
    void delete_nodes();
};

#endif // MY_LATTICE_BOLTZMANN_H
