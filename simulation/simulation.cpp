#include "simulation.h"

/**
 * Calculate the macro values with a force attached to it.
 * @todo macro integration force relation
 * @param a
 * @param previous_force
 * @return
 */
inline std::tuple<double, double, double> forcedSimulation::calculate_macro(array_t *a,
                                                                            array_t* previous_force,
                                                                            vector_t force) {
    // calculate rho ux and uy
    // shorthands
    int o = offset_node;
    double dt = parameters.dt;
    auto p = a->begin() + o;
    auto f = previous_force->begin();
    // rho, ux, uy
    // base terms sum
    /*
        ___      dt   ___
        ╲   f  + ── ⋅ ╲   F  = ρ
        ╱    i    2   ╱    i
        ‾‾‾           ‾‾‾
        i             i
     */
    double prefactor = dt/2;
    double rho = (p + 0).operator*() +
                 (p + 1).operator*() +
                 (p + 2).operator*() +
                 (p + 3).operator*() +
                 (p + 4).operator*() +
                 (p + 5).operator*() +
                 (p + 6).operator*() +
                 (p + 7).operator*() +
                 (p + 8).operator*();
    /*
    double f_term = weights(0) * (f+0).operator*() +
                    weights(1) * (f+1).operator*() +
                    weights(2) * (f+2).operator*() +
                    weights(3) * (f+3).operator*() +
                    weights(4) * (f+4).operator*() +
                    weights(5) * (f+5).operator*() +
                    weights(6) * (f+6).operator*() +
                    weights(7) * (f+7).operator*() +
                    weights(8) * (f+8).operator*();
    */
    double f_term = 0;
    rho += f_term*prefactor;
    // ux and uy
    /* Equation for the calculation of ux and uy
        ___           dt   ___
        ╲   f  ⋅ c  + ── ⋅ ╲   F  ⋅ c   = ρ ⋅ u
        ╱    i    i    2   ╱    i    ia
        ‾‾‾                ‾‾‾
         i                  i
     */
    //
    /*
    vector_t force_add = {0,0};
    for(int i = 0; i < CHANNELS; ++i) {
        vector_t velocity_channel = velocity_set.col(i);
        force_add += velocity_channel*(f+i).operator*()*weights(i);
    }
    force_add *= prefactor;
     */
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
    // ux and uy through
    // ux = (ux+force_add.x())/rho;
    // uy = (uy+force_add.y())/rho;
    ux = (ux + force.x())/rho;
    uy = (uy + force.y())/rho;
    return {rho, ux, uy};
}

/**
 * Streaming stepfunction .
 * @param a
 * @param list
 */
inline void forcedSimulation::streaming(array_t *a, std::vector<link_pointer> *list) {
    // just the sim
    for(int i = 1; i < CHANNELS; ++i) {
        // pointer magic
        auto origin = a->begin() + offset_node + i;
        (list->operator[](i-1) + offset_sim).operator*() = origin.operator*();
    }
}

/**
 * Collision term, the forcing is done in a subsequent step.
 * @param a
 * @param rho
 * @param ux
 * @param uy
 */
inline void forcedSimulation::collision(array_t *a, double rho, double ux, double uy) {
    // undrosed collision term
    int o = offset_node;
    double relaxation = 1/parameters.relaxation;
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
 * Forcing terms to calculate and appliy the force.
 * @todo alpha vs channel force relationship
 * @param n
 * @param write_to
 * @param ux
 * @param uy
 */
inline void forcedSimulation::forcing_terms(long pos, double ux, double uy) {
    // set some shorthands
    int o = offset_node;
    auto n = nodes[pos];
    auto write_to = forces[pos];
    auto force_a = force_alpha[pos];
    auto p = n->populations.begin() + o;
    rot_force->calculate_F_rotation(ux,uy,&n->position);
    //rot_force->calculate_F_circle(&n->position,0.0035,ux,uy);
    rot_force->calculate_F_i();
    force_a = rot_force->return_force_alpha();
    double prefactor = 1 - 1/(2*parameters.relaxation);
    for(int i = 0; i < CHANNELS; ++i) {
        double shorthand = rot_force->force_channels[i] * weights(i);
        // set into the force array
        write_to->operator[](i) = shorthand;
        // perform forcing step
        (p + i).operator*() += prefactor*shorthand;
    }
}

/**
 * Computes and spreads the lagrangian force over the channels.
 */
void forcedSimulation::compute_spread_lagrangian_force() {
    for(auto lm : lagrangian_markers) {
        // compute force
        vector_t force = compute_lagrangian_force();
        // spread force
        std::vector<handle_t> affected_nodes = all_nodes.ranging_key_translation(*lm,parameters.ibm_range);
        for(auto h : affected_nodes) {
            // influence the force
            h -= 1;
            double distance = (nodes[h]->position - *lm).norm();
            force_alpha[h] = kernel_3(distance,parameters.lattice_length) * force;
        }
    }
}

/**
 * Interpolates the lagrangian node position and forwards it in time.
 */
void forcedSimulation::interpolate_forward_lagrangian_force() {
    for(auto lm : lagrangian_markers) {
        vector_t u;
        std::vector<handle_t> affected_nodes = all_nodes.ranging_key_translation(*lm,parameters.ibm_range);
        for(auto h: affected_nodes) {
            vector_t expt = {1.201,12};
            h -= 1;
            double distance = (nodes[h]->position - *lm).norm();
            u = kernel_3(distance, parameters.lattice_length) * expt;
        }
    }
}
/**
 * Constructor.
 * @param c
 * @param g
 * @param f
 */
forcedSimulation::forcedSimulation(boundaryPointConstructor *c, nodeGenerator *g, goaForce *f) {
    boundary_points = c;
    node_generator = g;
    rot_force = f;
}

/**
 * Deconstructor.
 */
forcedSimulation::~forcedSimulation() {
    delete_nodes();
}

/**
 * Set the simulation parameters.
 * @param t
 */
void forcedSimulation::set_simulation_parameters(simulation_parameters_t t) {
    parameters = t;
    // fix the omega parameter to th new one
    parameters.relaxation = parameters.relaxation + parameters.dt/2;
}

/**
 * Inits the simulation class and allocates the nodes.
 */
void forcedSimulation::init() {
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
        for(auto link : node_info->links) {
            handle_t partner_handle = link.handle;
            int channel = link.channel;
            long array_position = long(partner_handle) - 1;
            auto link_p = nodes[array_position]->populations.begin() + channel;
            node->neighbors.push_back(link_p);
        }
    }
    for(auto n: nodes) {
        auto a = new array_t;
        a->resize(CHANNELS);
        a->setZero();
        forces.push_back(a);
        force_alpha.push_back({0,0});
    }
}

/**
 * Run function combines all the steps.
 * @param current_step
 */
void forcedSimulation::run(int current_step) {
    offset_sim = ((current_step +1) & 0x1) * 9;
    offset_node = (current_step & 0x1) * 9;
    for(long i = 0; i < nodes.size(); ++i) {
        // shorthands
        // node structure
        auto node = nodes[i];
        // force array
        auto force = forces[i];
        auto f_a = force_alpha[i];
        // use old force to update macro values
        auto [rho,ux, uy] = calculate_macro(&node->populations,force,f_a);
        // perform equilibrium step
        collision(&node->populations,rho,ux,uy);
        // add and calculate new force term
        forcing_terms(i,ux,uy);
        // streaming
        streaming(&node->populations,&node->neighbors);
    }
}

/**
 * Write out the node data into something python can plot.
 * @param write_to_file
 */
void forcedSimulation::get_data(bool write_to_file) {
    /// flowfields
    flowfield_t ux;
    flowfield_t uy;
    flowfield_t rho;
    long size_x = long(round(boundary_points->size.x()));
    long size_y = long(round(boundary_points->size.y()));
    ux.resize(size_x,size_y);
    uy.resize(size_x,size_y);
    rho.resize(size_x,size_y);
    // make sure to get the correct one
    for(int i = 0; i < nodes.size(); ++i) {
        auto node = nodes[i];
        auto force = forces[i];
        auto f_a = force_alpha[i];
        auto [rho_local,ux_local, uy_local] = calculate_macro(&node->populations,force,f_a);
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
 * Deletes and clears the vector data fields.
 */
void forcedSimulation::delete_nodes() {
    for (auto n : nodes) {
        delete n;
    }
    for (auto a : forces) {
        delete a;
    }
    nodes.clear();
    forces.clear();
}

/**
 * Hook to the inline method.
 * @param a
 * @param f
 * @return
 */
std::tuple<double,double,double> forcedSimulation::test_calcualte_macro(array_t *a, array_t *f) {
    vector_t t = {1,1};
    return calculate_macro(a,f,t);
}