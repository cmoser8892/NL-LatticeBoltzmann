#include "one_step_simulation.h"
#include "functions.h"
#include "helper_functions.h"
#include <iostream>

/// one step run class
/**
 * Function to calculate the rho, ux and uy values (macro values in the equilibrium function).
 * @param a
 * @return
 */
inline std::tuple<double, double, double> optimizedSimulation::calculate_macro(array_t *a) {
    // return als a struct
    // calculate rho ux and uy
    int o = offset_node;
    auto p = a->begin() + o;
    double rho = (p + 0).operator*() +
                 (p + 1).operator*() +
                 (p + 2).operator*() +
                 (p + 3).operator*() +
                 (p + 4).operator*() +
                 (p + 5).operator*() +
                 (p + 6).operator*() +
                 (p + 7).operator*() +
                 (p + 8).operator*();
    // ux + uy
    double ux = (((p + 1).operator*() +
                  (p + 5).operator*() +
                  (p + 8).operator*())-
                 ((p + 3).operator*() +
                  (p + 6).operator*()+
                  (p + 7).operator*()));
    double uy = (((p + 2).operator*() +
                  (p + 5).operator*() +
                  (p + 6).operator*())-
                 ((p + 4).operator*() +
                  (p + 7).operator*()+
                  (p + 8).operator*()));
    ux /= rho;
    uy /= rho;
    // return all the values
    return {rho, ux, uy};
}

/**
 * Calculates the forcing term.
 * @param n
 * @param ux
 * @param uy
 */
inline void optimizedSimulation::forcing_terms(oNode* n,double ux, double uy) {
    // set some shorthands
    int o = offset_node;
    auto p = n->populations.begin() + o;
    // precalculate the force
    //  rot_force->calculate_F_circle(&n->position,0.0035,ux,uy);
    rot_force->calculate_F_rotation(ux,uy,&n->position);
    // rot_force->calculate_F_circle(&n->position,0.0035);
    rot_force->calculate_F_i();
    for(int i = 0; i < CHANNELS; ++i) {
        (p + i).operator*() += rot_force->force_channels[i] * weights(i);
    }
}

/**
 * Simple bb based on the array only.
 * @param a
 */
inline void optimizedSimulation::bounce_back_moving(array_t *a) {
    // bb
    auto pointer = a->begin() + offset_sim;
    (pointer + 7).operator*() += -1.0/6 * parameters.u_wall;
    (pointer + 8).operator*() +=  1.0/6 * parameters.u_wall;
}

/**
 * Just array and linked list based streaming.
 * @param a
 * @param list
 */
inline void optimizedSimulation::streaming(array_t *a, std::vector<link_pointer>* list) {
    // just the sim
    for(int i = 1; i < CHANNELS; ++i) {
        // pointer magic
        auto origin = a->begin() + offset_node + i;
        (list->operator[](i-1) + offset_sim).operator*() = origin.operator*();
    }
}

/**
 * Optimized collision step.
 * @param a
 * @param rho
 * @param ux
 * @param uy
 */
inline void optimizedSimulation::collision(array_t* a, double rho, double ux, double uy) {
    // undrosed collision term
    int o = offset_node;
    double relaxation = parameters.relaxation;
    auto p = a->begin() + o;
    (p + 0).operator*() -= relaxation * ((p + 0).operator*() - weights.col(0).x()*rho*(1- 1.5*(ux*ux +uy*uy)));
    (p + 1).operator*() -= relaxation * ((p + 1).operator*() - weights.col(1).x()*rho*(1+ 3*ux+ 4.5*ux*ux- 1.5*(ux*ux +uy*uy)));
    (p + 2).operator*() -= relaxation * ((p + 2).operator*() - weights.col(2).x()*rho*(1+ 3*uy+ 4.5*uy*uy- 1.5*(ux*ux +uy*uy)));
    (p + 3).operator*() -= relaxation * ((p + 3).operator*() - weights.col(3).x()*rho*(1- 3*ux+ 4.5*ux*ux- 1.5*(ux*ux +uy*uy)));
    (p + 4).operator*() -= relaxation * ((p + 4).operator*() - weights.col(4).x()*rho*(1- 3*uy+ 4.5*uy*uy- 1.5*(ux*ux +uy*uy)));
    (p + 5).operator*() -= relaxation * ((p + 5).operator*() - weights.col(5).x()*rho*(1+ 3*ux+ 3*uy+ 9*ux*uy+ 3*(ux*ux +uy*uy)));
    (p + 6).operator*() -= relaxation * ((p + 6).operator*() - weights.col(6).x()*rho*(1- 3*ux+ 3*uy- 9*ux*uy+ 3*(ux*ux +uy*uy)));
    (p + 7).operator*() -= relaxation * ((p + 7).operator*() - weights.col(7).x()*rho*(1- 3*ux- 3*uy+ 9*ux*uy+ 3*(ux*ux +uy*uy)));
    (p + 8).operator*() -= relaxation * ((p + 8).operator*() - weights.col(8).x()*rho*(1+ 3*ux- 3*uy- 9*ux*uy+ 3*(ux*ux +uy*uy)));
}

/**
 * Constructor sets the sub-classes.
 * @param c
 * @param g
 */
optimizedSimulation::optimizedSimulation(boundaryPointConstructor *c, nodeGenerator *g) {
    boundary_points = c;
    node_generator = g;
    // only do if not given 0s
    if(g != nullptr && c != nullptr) {
        force = new gladrowForce(0.007,c->size);
    }
    else {
        std::cout << "You gave me nullptr, i will not work and crash" << std::endl;
    }

}

/**
 * Constructor with a force.
 * @param c
 * @param g
 * @param f
 */
optimizedSimulation::optimizedSimulation(boundaryPointConstructor *c, nodeGenerator *g, goaForce * f) {
    boundary_points = c;
    node_generator = g;
    rot_force = f;
}

/**
 * Deletes nodes.
 */
optimizedSimulation::~optimizedSimulation() {
    delete_nodes();
    delete force;
}

/**
 * Sets up the simulation parameters (given as a struct).
 * @param t
 */
void optimizedSimulation::set_simulation_parameters(simulation_parameters_t t) {
    parameters = t;
}

/**
 * Stream nodes out and performs the streaming step.
 * @param node
 */
void optimizedSimulation::streaming(oNode *node) {
    // loop through
    streaming(&node->populations,&node->neighbors);
}

/**
 * Performs the bounce back step on the top layer.
 * @param n
 */
void optimizedSimulation::bounce_back_moving(oNode *n) {
    if(n->boundary_type== BOUNCE_BACK_MOVING) {
        auto pointer = n->populations.begin() + offset_sim;
        (pointer + 7).operator*() += -1.0/6 * parameters.u_wall;
        (pointer + 8).operator*() +=  1.0/6 * parameters.u_wall;
    }
}

/**
 * One stop macro forcing term.
 * @param node
 */
void optimizedSimulation::one_step_macro_collision_forcing(oNode *node) {
    double relaxation = parameters.relaxation;
    int o = offset_node;
    // macro calc
    auto p = node->populations.begin() + o;
    // macro part
    // rho ux and uy inlined
    auto [rho,ux,uy] = calculate_macro(&node->populations);
    // collision
    (p + 0).operator*() -= relaxation *
                               ((p + 0).operator*() - weights.col(0).x()*rho*(1- 1.5*(ux*ux +uy*uy)))
                           + (0);
    (p + 1).operator*() -= relaxation *
                               ((p + 1).operator*() - weights.col(1).x()*rho*(1+ 3*ux+ 4.5*ux*ux- 1.5*(ux*ux +uy*uy)))
                           + (1.0/12*force->return_current_next_x(&node->position,1)) ;
    (p + 2).operator*() -= relaxation *
                               ((p + 2).operator*() - weights.col(2).x()*rho*(1+ 3*uy+ 4.5*uy*uy- 1.5*(ux*ux +uy*uy)))
                           + (1.0/12*force->return_current_next_y(&node->position,2));
    (p + 3).operator*() -= relaxation *
                               ((p + 3).operator*() - weights.col(3).x()*rho*(1- 3*ux+ 4.5*ux*ux- 1.5*(ux*ux +uy*uy)))
                           + (-1.0/12*force->return_current_next_x(&node->position,3));
    (p + 4).operator*() -= relaxation *
                               ((p + 4).operator*() - weights.col(4).x()*rho*(1- 3*uy+ 4.5*uy*uy- 1.5*(ux*ux +uy*uy)))
                           + (-1.0/12*force->return_current_next_y(&node->position,4));
    (p + 5).operator*() -= relaxation *
                               ((p + 5).operator*() - weights.col(5).x()*rho*(1+ 3*ux+ 3*uy+ 9*ux*uy+ 3*(ux*ux +uy*uy)))
                           + (1.0/12*force->return_current_next_x(&node->position,5)
                              + 1.0/12*force->return_current_next_y(&node->position,5));
    (p + 6).operator*() -= relaxation *
                               ((p + 6).operator*() - weights.col(6).x()*rho*(1- 3*ux+ 3*uy- 9*ux*uy+ 3*(ux*ux +uy*uy)))
                           + (-1.0/12*force->return_current_next_x(&node->position,6)
                              + 1.0/12*force->return_current_next_y(&node->position,6));
    (p + 7).operator*() -= relaxation *
                               ((p + 7).operator*() - weights.col(7).x()*rho*(1- 3*ux- 3*uy+ 9*ux*uy+ 3*(ux*ux +uy*uy)))
                           + (-1.0/12*force->return_current_next_x(&node->position,7) +
                              - 1.0/12*force->return_current_next_y(&node->position,7));
    (p + 8).operator*() -= relaxation *
                               ((p + 8).operator*() - weights.col(8).x()*rho*(1+ 3*ux- 3*uy- 9*ux*uy+ 3*(ux*ux +uy*uy)))
                           + (1.0/12*force->return_current_next_x(&node->position,8)
                              - 1.0/12*force->return_current_next_y(&node->position,8));
}

/**
 * One step for all the calculations necessary.
 * @param node
 * @param relaxation
 */
void optimizedSimulation::one_step_macro_collision(oNode* node, double relaxation) {
    // macro calc
    one_step_macro_collision(&node->populations);
}

/**
 * Performs the the macro and collision step in one go.
 * @param a
 */
void optimizedSimulation::one_step_macro_collision(array_t *a) {
    int o = offset_node;
    double relaxation = parameters.relaxation;
    // macro calc
    auto p = a->begin() + o;
    // rho ux and uy inlined
    auto [rho,ux,uy] = calculate_macro(a);
    // collision
    collision(a,rho,ux,uy);
}

/**
 * Inits the sim based on info in the node generator.
 */
void optimizedSimulation::init() {
    // set up nodes
    for(auto node_info : node_generator->node_infos) {
        auto n = new oNode(node_info->handle,velocity_set.cols(),node_info->boundary);
        // n->neighbors = node_info->links; // should copy everything not quite sure thou
        n->position = node_info->position;
        n->populations << equilibrium_2d(0,0,1) , equilibrium_2d(0,0,1);
        nodes.push_back(n);
    }
    // setup links
    for(int i = 0; i < node_generator->node_infos.size(); ++i) {
        // get them both in here
        auto node_info = node_generator->node_infos[i];
        auto node = nodes[i];
        //
        assert(node_info->links.size() > 8);
        for(auto link : node_info->links) {
            handle_t partner_handle = link.handle;
            int channel = link.channel;
            long array_position = long(partner_handle) - 1;
            auto link_p = nodes[array_position]->populations.begin() + channel;
            node->neighbors.push_back(link_p);
        }
    }
}

/**
 * Init the sub arrays too.
 */
void optimizedSimulation::init_sub_array() {
    // we push the information into the sub arrays
    for(auto n : nodes) {
        arrays_of_the_nodes.push_back(&n->populations);
        neighborhood_list.push_back(&n->neighbors);
        boundary.push_back(n->boundary_type);
    }
}

/**
 * Run the sim.
 * @param current_step
 */
void optimizedSimulation::run(int current_step ) {
    offset_sim = ((current_step +1) & 0x1) * 9;
    offset_node = (current_step & 0x1) * 9;
    for(auto n : nodes) {
        // macro and collision
        one_step_macro_collision(n,parameters.relaxation);
        // streaming
        streaming(n);
        // moving boundary
        bounce_back_moving(n);
    }
}

/**
 * Tests weather or not it makes sense to rip apart oNode (not really 3% improvement though)
 * @note Optimal if combined with a esoteric twist algorithm
 * @param current_step
 */
void optimizedSimulation::run_sub_array(int current_step) {
    offset_sim = ((current_step +1) & 0x1) * 9;
    offset_node = (current_step & 0x1) * 9;
    long range = nodes.size();
    for(int i = 0; i < range; ++i) {
        // shorthands
        auto population = arrays_of_the_nodes[i];
        auto pointer = neighborhood_list[i];
        auto bound = boundary[i];
        // functions
        one_step_macro_collision(population);
        streaming(population,pointer);
        // todo investigate cost of this statement / -> sort him out?!
        // todo prob one of the last possible optimizations left sort
        // so we dont have to load this var in memory
        if(bound == BOUNCE_BACK_MOVING) {bounce_back_moving(population);}
    }
}

/**
 * A forcing run.
 * @param current_step
 */
void optimizedSimulation::gladrow_force_run(int current_step) {
    offset_sim = ((current_step +1) & 0x1) * 9;
    offset_node = (current_step & 0x1) * 9;
    for(auto n : nodes) {
        // macro and collision
        //one_step_macro_collision(n,parameters.relaxation);
        one_step_macro_collision_forcing(n);
        // streaming
        streaming(n);
    }
}

/**
 * Implementation based on first order goa.
 * @param current_step
 */
void optimizedSimulation::forcing_run(int current_step) {
    offset_sim = ((current_step +1) & 0x1) * 9;
    offset_node = (current_step & 0x1) * 9;
    for(auto n : nodes) {
        // full inline now need to test with the constant circ force
        // macro
        auto [rho, ux, uy] = calculate_macro(&n->populations);
        // collision
        collision(&n->populations,rho,ux,uy);
        // force
        forcing_terms(n,ux,uy);
        // streaming
        streaming(&n->populations,&n->neighbors);
    }
}

/**
 * Writes data to file or prints it out.
 * @param write_to_file
 * @param orgiginalo
 */
void optimizedSimulation::get_data(bool write_to_file, point_t orgiginalo) {
    // we create a bunch of big arrays to store our data in
    // to be able to display like an array again
    flowfield_t ux;
    flowfield_t uy;
    flowfield_t rho;
    long size_x = long(round(orgiginalo.x()));
    long size_y = long(round(orgiginalo.y()));
    ux.resize(size_x,size_y);
    uy.resize(size_x,size_y);
    rho.resize(size_x,size_y);
    // make sure to get the correct one
    for(auto node: nodes) {
        auto [rho_local,ux_local, uy_local] = calculate_macro(&node->populations);
        ux(int(node->position(0)),int(node->position(1))) = ux_local;
        uy(int(node->position(0)),int(node->position(1))) = uy_local;
        rho(int(node->position(0)),int(node->position(1))) = rho_local;
    }
    // write to a file otherwise useless
    write_flowfield_data(&ux, "ux_data_file",write_to_file);
    write_flowfield_data(&uy, "uy_data_file",write_to_file);
    write_flowfield_data(&rho, "rho_data_file",write_to_file);
}

/**
 * Gets the data and writes it to file.
 * @param write_to_file
 */
void optimizedSimulation::get_data(bool write_to_file) {
    get_data(write_to_file,boundary_points->size);
}

/**
 * Frees up the memory again.
 */
void optimizedSimulation::delete_nodes() {
    for (auto n : nodes) {
        delete n;
    }
    nodes.clear();
}