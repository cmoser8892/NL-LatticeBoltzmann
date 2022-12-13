
#ifndef NL_LATTICEBOLTZMANN_NODE_GENERATOR_H
#define NL_LATTICEBOLTZMANN_NODE_GENERATOR_H

/**
 * class to generate/hold all the node points and connect them
 * second ingredient to the init
 * takes the boundary points as an ingredient
 */
#include "types.h"

typedef struct toLinks {
    int channel;
    handle_t handle;
}toLinks_t;

typedef struct nodePoint {
    handle_t handle;
    array_t position;
    nodeIdentifier_t type;
    std::vector<toLinks_t> links;
    boundaryType_t boundary; // ignored if node type is wet
}nodePoint_t;

class node_generator {
  private:
    // todo here
  public:
    std::vector<nodePoint_t*> node_infos;
};

#endif // NL_LATTICEBOLTZMANN_NODE_GENERATOR_H
