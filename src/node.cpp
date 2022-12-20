//
// Created by christoph on 13.12.22.
//
#include "node.h"

node::node(handle_t h, int dimensions, int channels, array_t pos) {
    handle = h;
    node_type = WET;
    rho = 1;
    data.setZero(channels);
    // todo make this variable ?!
    copy.setZero(channels);
    u.setZero(dimensions);
    position = pos;

}

boundary::boundary(handle_t h, int dimensions, int channels, array_t pos, boundaryType_t type)
    : node(h,dimensions,channels,pos){
    node_type = DRY;
    boundary_type = type;
}