#ifndef MY_LATTICE_BOLTZMANN_H
#define MY_LATTICE_BOLTZMANN_H

/*
 * Simulation classes:
 * The simulation classes are intentionally left at a their
 * final development stage to make it easier to understand the code as a whole.
 * Normally i would not do this as it produces more spaghetti code but I think it is reasonable to
 * do here.
 * New simulation classes are written when the current one can not without losing old features
 * be updated anymore or it just got to clunky as more and more ideas just got added on top.
 * Another benefit is backwards compatibility, which may be dropped at some point.
 * Everything not named simulation is considered to be out of data.
 * Heredity:
 * two_step_simulation (basic LB with NL works on just one array and a copy array)
 * one_step_simulation (optimized into one step algorithm
 *                      + most what is still doable without to much overhead)
 * simulation (forced second order integration)
 */

// includes
#include "boundary_point_generator.h"
#include "node.h"
#include "nodeGenerator.h"
#include "types.h"
#include "functions.h"
#include "forces.h"
// classen in gib
class forcedSimulation {
  private:
    simulation_parameters_t parameters;
    boundaryPointConstructor * boundary_points = nullptr;
    nodeGenerator* node_generator = nullptr;
    goaForce * rot_force = nullptr;
    // inlined core methods
    inline std::tuple<double, double, double> calculate_macro(array_t* a, array_t* previous_force);
    inline void streaming(array_t* a, std::vector<link_pointer> * list);
    inline void collision(array_t* a, double rho,double ux, double uy);
    inline void forcing_terms(oNode* n,array_t * w, double ux, double uy);
  public:
    // public main variables
    int offset_sim = 1;
    int offset_node = 0;
    // holding of the general nodes
    std::vector<oNode*> nodes;
    std::vector<array_t*> forces;
    // constructors
    forcedSimulation(boundaryPointConstructor* c,nodeGenerator* g, goaForce * f);
    ~forcedSimulation();
    // setters
    void set_simulation_parameters(simulation_parameters_t t);
    // run methods
    void init();
    void run(int current_step);
    // cleanup methods
    void get_data(bool write_to_file);
    void delete_nodes();
    // test functions
    std::tuple<double,double,double> test_calcualte_macro(array_t* a, array_t* f);
};
#endif // MY_LATTICE_BOLTZMANN_H
