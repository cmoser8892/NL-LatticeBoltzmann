#include "simulation.h"
#include "functions.h"
#include "helper_functions.h"
#include <iostream>

/// simulation run class
// private
void simulation::stream_links(node* n) {
    for(auto link : n->neighbors) {
        handle_t partner_handle = link->handle;
        int channel = link->channel;
        long array_position = long(partner_handle) - 1;
        // correct positioning prob
        nodes.at(array_position)->copy(channel)  = n->data(channel);
    }
}

void simulation::streaming_step_1() {
    for(auto node : nodes) {
        if(node->node_type == WET) {
            stream_links(node);
        }
    }
}

void simulation::streaming_step_2() {
    // basically just copy over data
    for(auto node : nodes) {
        node->data = node->copy;
    }
}

void simulation::bounce_back() {
    // aka a streaming step on boundary nodes only
    for(auto node : nodes) {
        if(node->node_type == DRY) {
            stream_links(node);
        }
    }
}

// public
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
        n->neighbors = node_info->links; // should copy everything not quite sure thou
        n->rho = 1;
        n->u.setZero();
        n->data = equilibrium(n);
        nodes.push_back(n);
    }
}

void simulation::run() {
    // run all substeps
    // moving wall missing i guess
    streaming_step_1();
    bounce_back();
    streaming_step_2();
    for(auto n : nodes) {
        macro(n);
    }
    for(auto n : nodes) {
        collision(n);
    }
}

void simulation::get_data() {
    flowfield_t ux;
    flowfield_t uy;
    flowfield_t rho;
    long size_x = long(round(boundary_points->size.x()));
    long size_y = long(round(boundary_points->size.y()));
    ux.resize(size_x,size_y);
    uy.resize(size_x,size_y);
    rho.resize(size_x,size_y);
    for(auto node: nodes) {
        // 2 methods that could be made into on, but for some indices
        write_ux(node,&ux);
        write_uy(node,&uy);
        write_rho(node,&rho);
    }
    // write to a file otherwise useless
    write_flowfield_data(&ux, "ux_flowfield");
    write_flowfield_data(&uy, "uy_flowfield");
    write_flowfield_data(&rho, "rho_flowfield");
}

std::vector<node*> simulation::access() {
    return nodes;
}


