#ifndef NL_LATTICEBOLTZMANN_ONE_STEP_SIMULATION_H
#define NL_LATTICEBOLTZMANN_ONE_STEP_SIMULATION_H

// includes
#include "boundary_point_generator.h"
#include "node.h"
#include "nodeGenerator.h"
#include "types.h"
#include "functions.h"
#include "forces.h"

/*
 * Notes:
 * The additional array for holding specific values of the sim bring
 * a performance increase of about 3% they were dropped in favor of a more
 * generalized structure.
 */
/**
 * Optimized version of the original basic Simulation class.
 * @note The additional array for holding specific values of the sim bring a performance increase of about 3% they were dropped in favor of a more generalized structure.
 */
class optimizedSimulation {
  private:
    simulation_parameters_t parameters; /**< simulation parameters, wall velocity and relaxation */
    boundaryPointConstructor * boundary_points = nullptr; /**< pointer to the previously constructed boundary point cloud */
    nodeGenerator* node_generator = nullptr; /**<  pointer to the to be used nodes */
    gladrowForce * force = nullptr; /**<  pointer to the rotation force calculation class */
    goaForce * rot_force = nullptr; /**<  pointer to the rotation force calculation class */
    // inlined core methods
    inline std::tuple<double, double, double> calculate_macro(array_t* a);
    inline void streaming(array_t* a, std::vector<link_pointer> * list);
    inline void collision(array_t* a, double rho,double ux, double uy);
    inline void forcing_terms(oNode* n, double ux, double uy);
    inline void bounce_back_moving(array_t* a);
  public:
    int offset_sim = 1; /**< offset 1 or 0 depending on the step */
    int offset_node = 0; /**<  offset 1 or 0 depending on the step (opposite to offset_sim) */
    // holding of the general nodes
    std::vector<oNode*> nodes; /**<  optimized nodes stored in a vector */
    // indirect pointers to the arrays
    std::vector<array_t*> arrays_of_the_nodes; /**<  pointers to the populations */
    // indirect pointers to the nl structure
    std::vector<std::vector<link_pointer>*> neighborhood_list; /**<  pointers to the neighborhood lists */
    std::vector<boundaryType_t> boundary; /**<  pointers to the type of boundaries */
    // methods
    optimizedSimulation(boundaryPointConstructor* c,nodeGenerator* g);
    optimizedSimulation(boundaryPointConstructor* c,nodeGenerator* g, goaForce * f);
    ~optimizedSimulation();
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
    void gladrow_force_run(int current_step);
    void forcing_run(int current_step);
    void get_data(bool write_to_file, point_t org);
    void get_data(bool write_to_file);
    void delete_nodes();
};

#endif // NL_LATTICEBOLTZMANN_ONE_STEP_SIMULATION_H
