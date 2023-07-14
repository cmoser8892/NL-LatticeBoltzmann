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

  public:
    // public main variables
    int offset_sim = 1; /**< offset 1 or 0 depending on the step */
    int offset_node = 0; /**<  offset 1 or 0 depending on the step (opposite to offset_sim) */
    std::vector<fNode*> nodes;
};

#endif // NL_LATTICEBOLTZMANN_LBM_SIMULATION_H
