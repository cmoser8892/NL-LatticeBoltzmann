#ifndef MY_LATTICE_BOLTZMANN_H
#define MY_LATTICE_BOLTZMANN_H

/*
 * Simulation classes:
 * The simulation classes are intentionally left at at their
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
/**
 * Currently worked on simulation class that supports forcing and is optimized for that.
 */
class forcedSimulation {
  private:
    simulation_parameters_t parameters; /**<  simulation parameters, wall velocity and relaxation */
    boundaryPointConstructor * boundary_points = nullptr; /**<  pointer to the previously constructed boundary point cloud */
    nodeGenerator* node_generator = nullptr; /**<  pointer to the to be used nodes */
    goaForce * rot_force = nullptr; /**<  pointer to the rotation force calculation class */
    markerIBM* original_markers;
    std::vector<point_t*> lagrangian_markers;
    rangingPointKeyHash all_nodes;
    coordinate_t size;
    // inlined core methods
    inline std::tuple<double, double, double> calculate_macro(array_t* a, array_t* previous_force,vector_t f);
    inline void streaming(array_t* a, std::vector<link_pointer> * list);
    inline void collision(array_t* a, double rho,double ux, double uy);
    inline void forcing_terms(long pos, double ux, double uy);
    inline void forcing_term_add(long pos, double ux, double uy, vector_t fa);
    void compute_spread_lagrangian_force();
    void interpolate_forward_lagrangian_force();
  public:
    // public main variables
    int offset_sim = 1; /**< offset 1 or 0 depending on the step */
    int offset_node = 0; /**<  offset 1 or 0 depending on the step (opposite to offset_sim) */
    // holding of the general nodes
    std::vector<oNode*> nodes; /**< optimized nodes stored in a vector */
    std::vector<array_t*> forces; /**<  force terms stored in a vector */
    std::vector<vector_t> force_alpha; /**< container of the force */
    // constructors
    forcedSimulation(boundaryPointConstructor* c,nodeGenerator* g, goaForce * f);
    forcedSimulation(nodeGenerator*g, goaForce*f, markerIBM* m, vector_t size);
    ~forcedSimulation();
    // setters
    void set_simulation_parameters(simulation_parameters_t t);
    // run methods
    void init();
    void run(int current_step);
    void run_ibm(int current_step);
    // cleanup methods
    void get_data(bool write_to_file);
    void delete_nodes();
    // test functions
    std::tuple<double,double,double> test_calcualte_macro(array_t* a, array_t* f);
};

#endif // MY_LATTICE_BOLTZMANN_H
