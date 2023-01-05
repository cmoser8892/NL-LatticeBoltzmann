#include "simulation.h"
#include "functions.h"
#include "helper_functions.h"
#include <iostream>

/// simulation run class
// private
int simulation::links_correct_channel(node * n, int link_channel) {
    int return_channel = -1;
    if(n->node_type == WET) {
        return_channel = link_channel;
    }
    else if( n-> node_type == DRY) {
        return_channel = switch_link_dimensions(link_channel);
    }
    else {
        throw std::invalid_argument("unknown node type");
    }
    return return_channel;
}

int simulation::switch_link_dimensions(int link_channel) {
    // aka hiding an ugly switch case
    int return_channel = -1;
    switch(link_channel) {
    case 1:
        return_channel = 3;
        break;
    case 2:
        return_channel = 4;
        break;
    case 3:
        return_channel = 1;
        break;
    case 4:
        return_channel = 2;
        break;
    case 5:
        return_channel = 7;
        break;
    case 6:
        return_channel = 8;
        break;
    case 7:
        return_channel = 5;
        break;
    case 8:
        return_channel = 6;
        break;
    default:
        throw std::invalid_argument("not possible dimension");
        break;
    }
    return return_channel;
}

void simulation::streaming_step_1() {
    for(auto node : nodes) {
        if(node->node_type == WET) {
            for(auto link : node->neighbors) {
                handle_t partner_handle = link->handle;
                int channel = link->channel;
                long array_position = long(partner_handle) - 1;
                // correct positioning prob
                nodes.at(array_position)->copy(channel)  = node->data(channel);
            }
        }
    }
}

void simulation::streaming_step_2() {
    // basically just copy over data
    for(auto node : nodes) {
        if(node->node_type == WET) {
            node->data = node->copy;
        }
    }
}

void simulation::bounce_back() {
    // aka a streaming step on boundary nodes only
    // when doing a bounce back it is crucical that all (wet) data is already in data and not in copy!!!
    for(auto node : nodes) {
        if(node->node_type == DRY) {
            for(auto link : node->neighbors) {
                // just for read ability will be optimized by the compiler
                handle_t partner_handle = link->handle;
                int link_channel = link->channel;
                int from_channel = links_correct_channel(node,link_channel);
                long array_position = long(partner_handle) - 1;
                // correct positioning
                double data = node->copy(from_channel);
                if(node->boundary_type == BOUNCE_BACK_MOVING) {
                    // apply the correct function to the channels
                    // no correction for different dry densities (rho_wall)
                    data += bb_switch_channel(from_channel,u_wall);
                }
                // directly write into the data
                nodes.at(array_position)->data(link_channel)  = data;
            }
        }
    }
}


// public
simulation::simulation(boundaryPointConstructor *c) {
    boundary_points = c;
}

simulation::simulation(boundaryPointConstructor *c, nodeGenerator *g) {
    boundary_points = c;
    node_generator = g;
}

void simulation::init() {
    // first initialize the node generator with the boundary points
    if(boundary_points == nullptr) {
        throw std::invalid_argument("no Boundary Points given");
    }
    /// todo dumb / setup for sliding lid
    double re = 1000; double base_length = boundary_points->size.x() - 2;
    u_wall = 0.1;
    relaxation = (2*re)/(6*base_length*u_wall+re);
    if(node_generator == nullptr) {
        // if the node generator hasnt run we have to run him
        node_generator = new nodeGenerator(boundary_points);
        node_generator->init();
    }
    // then rewrite the structure into the actual nodes
    for(auto node_info : node_generator->node_infos) {
        auto n = new node(node_info->handle,velocity_set.rows(),velocity_set.cols(),node_info->position,node_info->boundary);
        n->data.resize(velocity_set.cols());
        // todo make sure there are the right sizes and so on
        n->neighbors = node_info->links; // should copy everything not quite sure thou
        n->rho = 1;
        n->u.setZero();
        n->data = equilibrium(n);
        n->copy = n->data;
        nodes.push_back(n);
    }
}

void simulation::run() {
    // run all substeps
    // moving wall missing i guess
    streaming_step_1();
    streaming_step_2();
    bounce_back();
    for(auto n : nodes) {
        macro(n);
    }
    for(auto n : nodes) {
        collision(n,simulation::relaxation);
    }
}

void simulation::get_data(bool write_to_file) {
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
    write_flowfield_data(&ux, "ux_data_file",write_to_file);
    write_flowfield_data(&uy, "uy_data_file",write_to_file);
    write_flowfield_data(&rho, "rho_data_file",write_to_file);
}

std::vector<node*> simulation::access() {
    return nodes;
}


