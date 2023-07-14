#ifndef NL_LATTICEBOLTZMANN_LBM_SIMULATION_H
#define NL_LATTICEBOLTZMANN_LBM_SIMULATION_H

#include "node.h"
#include "nodeGenerator.h"
#include "types.h"
#include "functions.h"
#include "forces.h"

/**
 * Simulations with the Immersed boundary.
 */
class ibmSimulation {
  private:
    simulation_parameters_t parameters; /**<  simulation parameters, wall velocity and relaxation */
    nodeGenerator* node_generator = nullptr;
    rangingPointKeyHash markers_pkh;
    goaForce* rot_force = nullptr;
    // local access arrays
    inline std::tuple<double, double, double>  calcualte_macro(array_t* a, array_t* previous);
    inline void streaming(array_t* a, std::vector<link_pointer> *list);
    inline void collision(array_t* a, double rho);
    inline void forcing_term();
    inline void lbm_interpolate_forward_compute();
  public:
    // public main variables
    int offset_sim = 1; /**< offset 1 or 0 depending on the step */
    int offset_node = 0; /**<  offset 1 or 0 depending on the step (opposite to offset_sim) */
    std::vector<fNode*> nodes; /**< nodes-container */
    std::vector<marker*> marker; /**< marker-container */
    // constructor
    ibmSimulation(nodeGenerator* g, goaForce*f, markerIBM* m, vector_t size);
    ~ibmSimulation();
    // setters
    void set_simulation_parameters(simulation_parameters_t t);
    // runners
    void init();
    void run(int current_step);
    // past run methods
    void get_data(bool write_to_file);
    void delete_nodes();
    // test functions

};

#endif // NL_LATTICEBOLTZMANN_LBM_SIMULATION_H
