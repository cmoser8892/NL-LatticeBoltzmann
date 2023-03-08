#ifndef NL_LATTICEBOLTZMANN_NODE_H
#define NL_LATTICEBOLTZMANN_NODE_H

#include "types.h"

typedef struct nodePoint {
    handle_t handle;
    array_t position;
    nodeIdentifier_t type;
    std::vector<toLinks_t> links;
    boundaryType_t boundary; // ignored if node type is wet
}nodePoint_t;

class node {
  public:
    handle_t handle;
    nodeIdentifier_t node_type;
    boundaryType_t boundary_type;
    // data entries
    array_t* current_population;
    array_t*next_population;
    array_t population_even;
    array_t population_odd;
    std::vector<toLinks_t> neighbors;
    // macro values are local
    double rho;
    array_t u;
    array_t position;
    // pointer to functions no idea if they can work on the data and if i should
    // methods
    node(handle_t handle, int dimensions, int channels, array_t pos,boundaryType_t type);
};

class fusedNode {
    // handler info
    handle_t handle;
    boundaryType_t boundary_type;
    // data
    array_t* current_population;
    array_t* next_population;
    array_t population_even;
    array_t population_odd;
    // macro values
    double rho;
    array_t u;
    array_t position;
    fusedNode(handle_t handle, int dimensions, int channels, array_t pos,boundaryType_t type);
};

#endif // NL_LATTICEBOLTZMANN_NODE_H
