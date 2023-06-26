#include "node.h"

/**
 * Constructor sets all the variables.
 * @param h
 * @param dimensions
 * @param channels
 * @param pos
 * @param type
 */
node::node(handle_t h, int dimensions, int channels, array_t pos, boundaryType_t type ) {
    handle = h;
    boundary_type = type;
    if(type == NO_BOUNDARY) {
        node_type = WET;
    }
    else {
        node_type = DRY;
    }
    rho = 1;
    // 2 array approach
    population_even.setZero(channels);
    population_odd.setZero(channels);
    // set the pointers
    current_population = &population_even;
    next_population = &population_odd;
    u.setZero(dimensions);
    position = pos;
}

/**
 * oNode improvements.
 * @param h
 * @param channels
 * @param type
 */
oNode::oNode(handle_t h, int channels, boundaryType_t type) {
    handle = h;
    boundary_type = type;
    populations.setZero(2*channels);
}