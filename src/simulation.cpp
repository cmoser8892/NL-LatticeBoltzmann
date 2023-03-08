#include "simulation.h"
#include "functions.h"
#include "helper_functions.h"
#include <iostream>

/// simulation run class
// private

/// todo make the copy field non sticky aka delete stuff inbetween?!
/// todo no indirect access
/**
 * @fn void simulation::streaming_step_1()
 * @brief streaming is done via two steps first we put info in the copy field in the node
 */
void simulation::streaming_step_1() {
    // streaming is only for wet nodes!!!!
    for(auto node : nodes) {
        if(node->node_type == WET) {
            // very important so that we later not rewrite the value back to equilibrium
            for(auto link : node->neighbors) {
                handle_t partner_handle = link.handle;
                int channel = link.channel;
                long array_position = long(partner_handle) - 1;
                // correct positioning prob
                nodes.at(array_position)->population_odd(channel)  = node->population_even(channel);
            }
        }
    }
}

/**
 * @fn void simulation::streaming_step_2()
 * @brief we copy the data from copy over to the data field
 */
void simulation::streaming_step_2() {
    // basically just copy over data, importantly not channel 0!
    // streaming is only for wet nodes!!!!
    for(auto node : nodes) {
        if(node->node_type == WET) {
            for(int i = 1; i < CHANNELS; ++i) {
                node->population_even(i) = node->population_odd(i);
            }
        }
    }
}

/**
 * @fn void simulation::bounce_back()
 * @brief does the bounce back for all the dry nodes
 */
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
                double data = node->population_odd(from_channel);
                if(node->boundary_type == BOUNCE_BACK_MOVING) {
                    // apply the correct function to the channels
                    // no correction for different dry densities (rho_wall)
                    data += bb_switch_channel(from_channel,parameters.u_wall);
                }
                // directly write into the data
                nodes.at(array_position)->population_even(link_channel)  = data;
            }
        }
    }
}

/**
 * @fn void simulation::collisions()
 * @brief collision step
 */
void simulation::collisions() {
    // todo investigate singling out the equilibrium calculation
    double relax = parameters.relaxation;
    for(auto node: nodes) {
        /*
        // prob causes the cpu to stall/branch miss-prediction is more costly than not calculating
        // sometimes also doesnt do anything lol
         if(node->node_type == DRY)
            continue;
        */
        // calculate the macro values first
        macro(node);
        // convenience programming
        double ux = node->u(0);
        double uy = node->u(1);
        double rho = node->rho;
        // unroll the collision function also messy for a reason optimizes better
        node->population_even(0) -= relax * (node->population_even(0) - weights.col(0).x()*rho*(1- 1.5*(ux*ux +uy*uy)));
        node->population_even(1) -= relax * (node->population_even(1) - weights.col(1).x()*rho*(1+ 3*ux+ 4.5*ux*ux- 1.5*(ux*ux +uy*uy)));
        node->population_even(2) -= relax * (node->population_even(2) - weights.col(2).x()*rho*(1+ 3*uy+ 4.5*uy*uy- 1.5*(ux*ux +uy*uy)));
        node->population_even(3) -= relax * (node->population_even(3) - weights.col(3).x()*rho*(1- 3*ux+ 4.5*ux*ux- 1.5*(ux*ux +uy*uy)));
        node->population_even(4) -= relax * (node->population_even(4) - weights.col(4).x()*rho*(1- 3*uy+ 4.5*uy*uy- 1.5*(ux*ux +uy*uy)));
        node->population_even(5) -= relax * (node->population_even(5) - weights.col(5).x()*rho*(1+ 3*ux+ 3*uy+ 9*ux*uy+ 3*(ux*ux +uy*uy)));
        node->population_even(6) -= relax * (node->population_even(6) - weights.col(6).x()*rho*(1- 3*ux+ 3*uy- 9*ux*uy+ 3*(ux*ux +uy*uy)));
        node->population_even(7) -= relax * (node->population_even(7) - weights.col(7).x()*rho*(1- 3*ux- 3*uy+ 9*ux*uy+ 3*(ux*ux +uy*uy)));
        node->population_even(8) -= relax * (node->population_even(8) - weights.col(8).x()*rho*(1+ 3*ux- 3*uy- 9*ux*uy+ 3*(ux*ux +uy*uy)));
        /*
        // channel version
        for(int i = 0; i < CHANNELS; ++i) {
            double cx = velocity_set.col(i).x();
            double cy = velocity_set.col(i).y();
            double w = weights.col(i).x();
            node->data(i) -= relax * (node->data(i) - w*rho*(1 + 3*cx*ux + 3*cy*uy + 4.5*cx*cx*ux*ux + 9*cx*ux*cy*uy + 4.5*cy*cy*uy*uy - 1.5*ux*ux - 1.5*uy*uy));
        }
         */
        // og
        //node->data -= relax * (node->data - equilibrium(node));
    }
}
/**
 * @fn
 * @brief
 * @param node
 */
void simulation::fused_streaming(node *node) {
    // loop through
    for(int i = 1; i < CHANNELS; ++i) {
        auto link = node->neighbors.at(i-1);
        handle_t partner_handle = link.handle;
        int channel = link.channel;
        long array_position = long(partner_handle) - 1;
        // correct positioning prob
        nodes.at(array_position)->next_population->operator()(channel) = node->current_population->operator()(i);
    }
}

/**
 * @fn simulation::simulation(boundaryPointConstructor *c)
 * @brief constructor
 * @param c
 */
simulation::simulation(boundaryPointConstructor *c) {
    boundary_points = c;
}

/**
 * @fn simulation::simulation(boundaryPointConstructor *c, nodeGenerator *g)
 * @brief constructor
 * @param c
 * @param g
 */
simulation::simulation(boundaryPointConstructor *c, nodeGenerator *g) {
    boundary_points = c;
    node_generator = g;
}

/**
 * @fn simulation::~simulation()
 * @brief deconstructor, deletes all the nodes to prevent memory leaks
 */
simulation::~simulation() {
    delete_nodes();
}
void simulation::set_simulation_parameters(const simulation_parameters_t t) {
    parameters = t;
}

/**
 * @fn void simulation::init()
 * @brief initialize the nodes, without the calling this function once sim wont work
 */
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
    // then rewrite the structure into the actual nodes
    for(auto node_info : node_generator->node_infos) {
        auto n = new node(node_info->handle,velocity_set.rows(),velocity_set.cols(),node_info->position,node_info->boundary);
        n->population_even.resize(CHANNELS);
        // todo make sure there are the right sizes and so on
        n->neighbors = node_info->links; // should copy everything not quite sure thou
        n->rho = 1;
        n->u.setZero();
        n->population_even = equilibrium(n);
        n->population_odd = n->population_even;
        nodes.push_back(n);
    }
}

/**
 *  @fn void simulation::fused_init()
 *  @brief fused init
 */
void simulation::fused_init() {
    // first initialize the node generator with the boundary points
    if(boundary_points == nullptr) {
        throw std::invalid_argument("no Boundary Points given");
    }
    if(node_generator == nullptr) {
        throw std::invalid_argument("no Node-generator given");
    }
    // then rewrite the structure into the actual nodes
    for(auto node_info : node_generator->node_infos) {
        auto n = new node(node_info->handle,velocity_set.rows(),velocity_set.cols(),node_info->position,node_info->boundary);
        n->neighbors = node_info->links; // should copy everything not quite sure thou
        n->rho = 1;
        n->u.setZero();
        n->population_even = equilibrium(n);
        n->population_odd = n->population_even;
        nodes.push_back(n);
    }
}
/**
 * @fn void simulation::run()
 * @brief on sim run
 */
void simulation::run() {
    // todo reduce full runs through the nodes
    // run all substeps
    // moving wall missing i guess
    streaming_step_1();
    streaming_step_2();
    bounce_back();
    /*
    for(auto n: nodes) {
        macro(n);
    }
    for(auto n : nodes) {
        collision(n,parameters.relaxation);
    }
     */
    collisions(); // also does macro
}

/**
 * @fn void simulation::fused_run()
 * @brief this is a one step implementation of lb using the already described nodes
 */
void simulation::fused_run() {
    for(auto node : nodes) {
        // streaming and bb
        fused_streaming(node);
        // macro and collision
        fused_macro(node);
        fused_collision(node,parameters.relaxation);
        // switchero
        array_t * temp = node->current_population;
        node->current_population = node->next_population;
        node->next_population = temp;
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
    long size_x = long(round(original_size.x()));
    long size_y = long(round(original_size.y()));
    ux.resize(size_x,size_y);
    uy.resize(size_x,size_y);
    rho.resize(size_x,size_y);
    // make sure to get the correct
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

/**
 * @fn void simulation::delete_nodes()
 * @brief deletes the nodes as it says
 */
void simulation::delete_nodes() {
    for (auto n : nodes) {
        delete n;
    }
}