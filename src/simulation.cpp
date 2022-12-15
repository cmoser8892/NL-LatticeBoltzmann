#include "simulation.h"
#include "functions.h"
#include "helper_functions.h"
#include <iostream>

/// simulation run class
simulation::simulation(boundaryPointConstructor *c) {
    boundary_points = c;
}

void simulation::init() {
    // first initialize the node generator with the boundary points
    node_generator = new nodeGenerator(boundary_points);
    node_generator->init();
    // then rewrite the structure into the actual nodes
    for(auto node_info : node_generator->node_infos) {
        auto n = new node(node_info->handle,velocity_set.rows(),velocity_set.cols(),node_info->position);
        n->data.resize(velocity_set.cols());
        // todo make sure there are the right sizes and so on
        nodes.push_back(n);
    }
    // apply channel knowledge
    // todo
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


