//
// Created by christoph on 13.12.22.
//

#ifndef NL_LATTICEBOLTZMANN_NODE_H
#define NL_LATTICEBOLTZMANN_NODE_H

#include "types.h"

class node {
  public:
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
    node(int dimensions, int channels, array_t pos);
};

class boundary : public node {
    // a boundary is a special kind of node
    // they are stored
  public:
    boundaryType_t boundary_type;
    boundary(int dimensions, int channels, array_t pos, boundaryType_t type);
};

#endif // NL_LATTICEBOLTZMANN_NODE_H
