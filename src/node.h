
#ifndef NL_LATTICEBOLTZMANN_NODE_H
#define NL_LATTICEBOLTZMANN_NODE_H

#include "types.h"

class node {
  public:
    handle_t handle;
    nodeIdentifier_t node_type;
    array_t data;
    array_t copy;
    std::vector<node*> neighbors;
    // macro values are local
    double rho;
    array_t u;
    array_t position;
    // pointer to functions no idea if they can work on the data and if i should
    // methods
    node(handle_t handle, int dimensions, int channels, array_t pos);
};

class boundary : public node {
    // a boundary is a special kind of node
    // they are stored
  public:
    boundaryType_t boundary_type;
    boundary(handle_t h, int dimensions, int channels, array_t pos, boundaryType_t type);
};

#endif // NL_LATTICEBOLTZMANN_NODE_H
