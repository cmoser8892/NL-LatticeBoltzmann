#include "node.h"

/**
 * @fn node::node(handle_t h, int dimensions, int channels, array_t pos, boundaryType_t type )
 * @brief constructor sets all the variables
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
    u.setZero(dimensions);
    position = pos;
}