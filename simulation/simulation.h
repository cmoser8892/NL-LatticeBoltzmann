#ifndef MY_LATTICE_BOLTZMANN_H
#define MY_LATTICE_BOLTZMANN_H

// includes
#include "boundary_point_generator.h"
#include "node.h"
#include "nodeGenerator.h"
#include "types.h"
#include "functions.h"
#include "forces.h"

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
    // todo need to combine the different force at some point
    gladrowForce * force = nullptr;
    goaForce * rot_force = nullptr;
    // inlined core methods
    inline std::tuple<double, double, double> calculate_macro(array_t* a);
    inline void streaming(array_t* a, std::vector<link_pointer> * list);
    inline void collision(array_t* a, double rho,double ux, double uy);
    inline void forcing_terms(oNode* n, double ux, double uy);
    inline void bounce_back_moving(array_t* a);
  public:
    int offset_sim = 1;
    int offset_node = 0;
    // holding of the general nodes
    std::vector<oNode*> nodes;
    // indirect pointers to the arrays
    std::vector<array_t*> arrays_of_the_nodes;
    // indirect pointers to the nl structure
    std::vector<std::vector<link_pointer>*> neighborhood_list;
    std::vector<boundaryType_t> boundary;
    // methods
    oSimu(boundaryPointConstructor* c,nodeGenerator* g);
    oSimu(boundaryPointConstructor* c,nodeGenerator* g, goaForce * f);
    ~oSimu();
    void set_simulation_parameters(simulation_parameters_t t);
    // test methods
    void streaming(oNode* n);
    void bounce_back_moving(oNode* n);
    void one_step_macro_collision_forcing(oNode* n);
    void one_step_macro_collision(oNode* n, double relaxation);
    void one_step_macro_collision(array_t* a);
    // other stuff
    void init();
    void init_sub_array();
    void run(int current_step);
    void run_sub_array(int current_step);
    void current_run(int current_step);
    void forcing_run(int current_step);
    void get_data(bool write_to_file, point_t org);
    void get_data(bool write_to_file);
    void delete_nodes();
};

#endif // MY_LATTICE_BOLTZMANN_H
