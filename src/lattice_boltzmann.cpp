#include "lattice_boltzmann.h"

node::node(int dimensions, int channels, array_t position) {

    data.resize(channels);
    u.resize(dimensions);
    position.resize(dimensions);
    position = position;
    // neighbors are not setup here after creation
}
