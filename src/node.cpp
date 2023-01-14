//
// Created by christoph on 13.12.22.
//
#include "node.h"

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
    data.setZero(channels);
    // todo make this variable ?!
    copy.setZero(channels);
    u.setZero(dimensions);
    position = pos;
}