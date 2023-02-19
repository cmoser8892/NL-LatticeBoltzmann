#include "simulation.h"
#include "functions.h"
#include "helper_functions.h"
#include <iostream>

/// simulation run class
// private
void simulation::streaming_step_1() {
    for(auto node : nodes) {
        if(node->node_type == WET) {
            // very important so that we later not rewrite the value back to equilibrium
            for(auto link : node->neighbors) {
                handle_t partner_handle = link.handle;
                int channel = link.channel;
                long array_position = long(partner_handle) - 1;
                // correct positioning prob
                nodes.at(array_position)->copy(channel)  = node->data(channel);
            }
        }
    }
}

void simulation::streaming_step_2() {
    // basically just copy over data, importantly not channel 0!
    for(auto node : nodes) {
        if(node->node_type == WET) {
            for(int i = 1; i < node->data.size(); ++i) {
                node->data(i) = node->copy(i);
            }
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
                handle_t partner_handle = link.handle;
                int link_channel = link.channel;
                int from_channel = switch_link_dimensions(link_channel);
                long array_position = long(partner_handle) - 1;
                // correct positioning
                double data = node->copy(from_channel);
                if(node->boundary_type == BOUNCE_BACK_MOVING) {
                    // apply the correct function to the channels
                    // no correction for different dry densities (rho_wall)
                    data += bb_switch_channel(from_channel,parameters.u_wall);
                }
                // directly write into the data
                nodes.at(array_position)->data(link_channel)  = data;
            }
        }
    }
}

void simulation::collisions() {
    double relax = parameters.relaxation;
    for(auto node: nodes) {
        node->data -= relax * (node->data - table->equilibrium(node));
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

simulation::~simulation() {
    for (auto n : nodes) {
        delete n;
    }
}
void simulation::set_simulation_parameters(const simulation_parameters_t t) {
    parameters = t;
}


void simulation::init() {
    // first initialize the node generator with the boundary points
    if(boundary_points == nullptr) {
        throw std::invalid_argument("no Boundary Points given");
    }
    if(node_generator == nullptr) {
        // if the node generator hasnt run we have to run him
        node_generator = new nodeGenerator(boundary_points);
        node_generator->init();
    }
    if(table == nullptr) {
        // setup lookup
        table = new lookup(parameters.lookup_bits,
                           parameters.lookup_floor,
                           parameters.lookup_ceiling,
                           false);
        table->set_bypass(parameters.bypass_lookup);
    }
    // then rewrite the structure into the actual nodes
    for(auto node_info : node_generator->node_infos) {
        auto n = new node(node_info->handle,velocity_set.rows(),velocity_set.cols(),node_info->position,node_info->boundary);
        n->data.resize(velocity_set.cols());
        // todo make sure there are the right sizes and so on
        n->neighbors = node_info->links; // should copy everything not quite sure thou
        n->rho = 1;
        n->u.setZero();
        n->data = table->equilibrium(n);
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
    /*
    for(auto n : nodes) {
        collision(n,parameters.relaxation);
    } */
    collisions();
    for(auto n : nodes) {
        //debug_node(n, true);
    }
}
/**
 * @fn void simulation::get_data(bool write_to_file)
 * @brief puts the ux, uy and rho of individual nodes into a flow-field for visualization
 * @param write_to_file write to file or to cout
 */
void simulation::get_data(bool write_to_file, point_t original_size) {
    flowfield_t ux;
    flowfield_t uy;
    flowfield_t rho;
    std::cout << "Recalcs " <<table->recalc << std::endl;
    std::cout << "Table-hits " << table->table_hits << std::endl;
    std::cout << "Non-hits " << table->non_hits << std::endl;
    long size_x = long(round(original_size.x()));
    long size_y = long(round(original_size.y()));
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


