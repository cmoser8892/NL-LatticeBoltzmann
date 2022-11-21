#include "lattice_boltzmann.h"

/// nodes
// constructor should be the only thing needed
node::node(int dimensions, int channels, array_t position) {
    data.resize(channels);
    copy.resize(channels);
    u.resize(dimensions);
    position.resize(dimensions);
    position = position;
    // neighbors are not setup here after creation
    equilibrium_func = nullptr;
    streaming_func = nullptr;
    collision_func = nullptr;
    macro_func = nullptr;
}

/// simulation run class
void simulation::init() {

}

void simulation::run() {

}



