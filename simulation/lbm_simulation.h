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
    coordinate_t size;
    nodeGenerator* node_generator = nullptr;
    goaForce* rot_force = nullptr;
    markerIBM* original_markers = nullptr;
    rangingPointKeyHash markers_pkh;
    // local access arrays
    inline std::tuple<double, double, double> calculate_macro(array_t * a);
    inline std::tuple<double, double, double>  calculate_macro_force(array_t* a, array_t* previous);
    inline void streaming(array_t* a, std::vector<link_pointer> *list);
    inline void collision(array_t* a, double rho, double ux, double uy);
    vector_t calculate_rotation_force(point_t* pos, vector_t* velocity);
    inline void forcing_term(fNode* n, vector_t* force);
    vector_t aggregate_force(std::vector<handle_t> *handles, point_t* pos);
    void distribute_velocity(std::vector<handle_t> *handles, point_t* pos, vector_t * v);
    void propagate_calculate_force_marker();
  public:
    // public main variables
    int offset_sim = 1; /**< offset 1 or 0 depending on the step */
    int offset_node = 0; /**<  offset 1 or 0 depending on the step (opposite to offset_sim) */
    std::vector<fNode*> nodes; /**< nodes-container */
    std::vector<marker*> markers; /**< marker-container */
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
    void delete_containers();
    // test functions

};

#endif // NL_LATTICEBOLTZMANN_LBM_SIMULATION_H
