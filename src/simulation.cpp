#include "simulation.h"
#include "functions.h"
#include "helper_functions.h"
#include <iostream>

/// nodes

/// simulation run class
void simulation::init(int six, int siy) {

}

void simulation::run() {

}

void simulation::get_data() {
    flowfield_t ux;
    flowfield_t uy;
    flowfield_t rho;
    ux.resize(size_x,size_y);
    uy.resize(size_x,size_y);
    rho.resize(size_x,size_y);
    for(auto node: nodes) {
        // 2 methods that could be made into on, but for some indices
        write_ux(node,&ux);
        write_uy(node,&uy);
    }
    // write to a file otherwise useless
    write_flowfield_data(&ux, "ux_flowfield");
    write_flowfield_data(&uy, "uy_flowfield");
    write_flowfield_data(&rho, "rho_flowfield");
}

std::vector<node*> simulation::access() {
    return nodes;
}


