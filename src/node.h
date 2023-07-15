#ifndef NL_LATTICEBOLTZMANN_NODE_H
#define NL_LATTICEBOLTZMANN_NODE_H

#include "types.h"

/**
 * Lightweight node used for construction without all the data fields.
 */
typedef struct nodePoint {
    handle_t handle; /**<  Handle of the node */
    array_t position; /**<  Position */
    nodeIdentifier_t type; /**<  Dry or wet */
    std::vector<toLinks_t> links; /**<  Neighbors of that node */
    boundaryType_t boundary; /**<  Ignored if node type is wet */
}nodePoint_t;

// node class mostly replaced with oNode in actual code, still used for testing
/**
 * Node class (unoptimized).
 */
class node {
  public:
    handle_t handle; /**< Handle to the node */
    nodeIdentifier_t node_type; /**<  Type of node wet or dry */
    boundaryType_t boundary_type; /**<  Boundary type */
    // data entries
    array_t* current_population; /**<  current population pointer */
    array_t* next_population; /**<  next population pointer */
    array_t population_even; /**<  even step population */
    array_t population_odd; /**<  odd step population */
    std::vector<toLinks_t> neighbors; /**<  Node neighbors */
    // macro values are local
    double rho; /**<  rho */
    array_t u; /**<  velocity */
    array_t position; /**<  Position of the node */
    // pointer to functions no idea if they can work on the data and if i should
    // methods
    node(handle_t handle, int dimensions, int channels, array_t pos,boundaryType_t type);
};

/**
 * Optimized node with less data and therefore better handling.
 */
class oNode {
  public:
    // handle and boundary
    handle_t handle; /**< Handle to the node */
    point_t position; /**< Position of the node */
    vector_t velocity; /**< Velocity */
    boundaryType_t boundary_type; /**< Type of boundary (Bounce back, Moving, etc.) */
    // 2xChannels
    array_t populations; /**< 2x9 population */
    std::vector<link_pointer> neighbors; /**< The neighbors of the node */
    // constructor
    oNode(handle_t, int channels, boundaryType_t type);
};

/**
 * fNode for ibm
 */
class fNode {
  public:
    handle_t handle;              /**< Handle to the node */
    point_t position;             /**< Position of the node */
    vector_t velocity;            /**< Velocity */
    boundaryType_t boundary_type; /**< Type of boundary (Bounce back, Moving, etc.) */
    //
    array_t populations;                 /**< 2x9 population */
    array_t forces;
    std::vector<link_pointer> neighbors; /**< The neighbors of the node */
    // constructor
    fNode(handle_t, int channels, boundaryType_t type);
};

/**
 * ibm marker
 */
class marker {
  public:
    point_t position; /**< Location of the marker*/
    vector_t force; /**< Force affecting the marker */
    point_t original_position; /**< Original location of the marker */
    explicit marker(point_t pos);
};

#endif // NL_LATTICEBOLTZMANN_NODE_H
