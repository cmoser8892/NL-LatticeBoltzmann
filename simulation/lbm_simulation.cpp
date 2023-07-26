# include "lbm_simulation.h"

std::tuple<double, double, double> ibmSimulation::calculate_macro(array_t* a) {
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
    ux = ux/rho;
    uy = uy/rho;
    return {rho, ux, uy};
}

std::tuple<double, double, double> ibmSimulation::calculate_macro_force(array_t *a, array_t *previous_force) {
    // calculate rho ux and uy
    // shorthands
    int o = offset_node;
    double dt = parameters.dt;
    auto p = a->begin() + o;
    auto f = previous_force->begin();
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
    double f_term = weights(0) * (f+0).operator*() +
                    weights(1) * (f+1).operator*() +
                    weights(2) * (f+2).operator*() +
                    weights(3) * (f+3).operator*() +
                    weights(4) * (f+4).operator*() +
                    weights(5) * (f+5).operator*() +
                    weights(6) * (f+6).operator*() +
                    weights(7) * (f+7).operator*() +
                    weights(8) * (f+8).operator*();
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
    vector_t force_add = {0,0};
    for(int i = 0; i < CHANNELS; ++i) {
        vector_t velocity_channel = velocity_set.col(i);
        force_add += velocity_channel*(f+i).operator*()*weights(i);
    }
    force_add *= prefactor;
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
    ux = (ux+force_add.x())/rho;
    uy = (uy+force_add.y())/rho;
    return {rho, ux, uy};
}
void ibmSimulation::streaming(array_t *a, std::vector<link_pointer> *list) {
    // just the sim
    for(int i = 1; i < CHANNELS; ++i) {
        // pointer magic
        auto origin = a->begin() + offset_node + i;
        (list->operator[](i-1) + offset_sim).operator*() = origin.operator*();
    }
}
void ibmSimulation::collision(array_t *a, double rho, double ux, double uy) {
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

vector_t ibmSimulation::calculate_rotation_force(point_t* pos, vector_t *velocity) {
    // rot_force->calculate_F_rotation(velocity->x(),velocity->y(),pos);
    rot_force->calculate_F_circle(pos,0.0001,velocity->x(),velocity->y());
    return rot_force->return_force_alpha();
}

void ibmSimulation::forcing_term(fNode* n,vector_t* force) {
    // calculate the forcing term and use add the aggregate boundary force
    // calculate fi and correct the original collision term
    // array_t force_i = calculate_truncation_array()
    int o = offset_node;
    auto p = n->populations.begin() + o;
    double prefactor = 1 - 1/(2*parameters.relaxation);
    // calculate the terms
    array_t f = calculate_truncation_array(force,&n->velocity);
    // apply to the post collision terms
    for(int i = 0; i < CHANNELS; ++i) {
        double shorthand = f[i] * weights(i);
        // put it in the save array
        n->forces[i] = shorthand;
        // perform forcing step
        (p + i).operator*() += prefactor*shorthand;
    }
}

vector_t ibmSimulation::aggregate_force(std::vector<handle_t> *handles, point_t* pos) {
    vector_t returns ={0,0};
    // we look for markers in the vicinity of our node and add up the total force
    for(auto h : *handles) {
        h-= 1;
        auto marker = markers[h];
        vector_t r = (*pos - marker->position);
        returns += (*kernel_function)(&r)* marker->force;
    }
    return returns;
}

void ibmSimulation::distribute_velocity(std::vector<handle_t> *handles,  point_t* pos, vector_t * v) {
    // we distribute and put the velocity back into the makers
    for(auto h : *handles) {
        h-= 1;
        auto marker = markers[h];
        vector_t r = (*pos - marker->position);
        marker->velocity += (*kernel_function)(&r) * (*v);
    }
}

void ibmSimulation::propagate_calculate_force_marker() {
    // we do this after dealing with the force markers
    for(auto m : markers) {
        // update the position
        m->position += m->velocity*parameters.dt;
        // calculate a new force based on the orginal position and so on
        vector_t delta_r = m->position - m->original_position;
        m->force = -parameters.k * (parameters.mean_marker_distance/parameters.lattice_length) * delta_r;
        // rot_force->calculate_F_rotation(m->velocity.x(),m->velocity.y(),&m->position);
        // m->force += rot_force->return_force_alpha();
        // set the velocity back to 0 so that we can add back up
        m->velocity.setZero();
    }
}

double ibmSimulation::kernel_function_call(point_t *p) {
    vector_t temp = *p/2;
    return (*kernel_function)(&temp);
}

ibmSimulation::ibmSimulation(nodeGenerator *g, goaForce *f, markerIBM *m, vector_t s) {
    node_generator = g;
    rot_force = f;
    original_markers = m;
    size.x = s.x();
    size.y = s.y();
}

ibmSimulation::~ibmSimulation() {
    delete_containers();
}
void ibmSimulation::set_simulation_parameters(simulation_parameters_t t) {
    parameters = t;
    // fix the omega parameter to th new one
    parameters.relaxation = parameters.relaxation + parameters.dt/2;
    // set the correct kernel function
    switch (t.kernel_in_use){
    case KERNEL_A: {
        kernel_function = &kernel_A_2d;
        break;
    }
    case KERNEL_B: {
        kernel_function = &kernel_B_2d;
        break;
    }
    case KERNEL_C: {
        kernel_function = &kernel_C_2d;
        break;
    }
    default:
        std::cerr << "Unknown kernel Type (IBM)" << std::endl;
        break;
    }
    // set the correct range
    // ATTENTION should not be set outside of this function
    parameters.ibm_range =kernel_id_to_lattice_search(t.kernel_in_use);
}
void ibmSimulation::init() {
    // setup the basic nodes
    for(auto node_info : node_generator->node_infos) {
        auto n = new fNode(node_info->handle,velocity_set.cols(),node_info->boundary);
        // n->neighbors = node_info->links; // should copy everything not quite sure thou
        n->position = node_info->position;
        n->populations << equilibrium_2d(0,0,1) , equilibrium_2d(0,0,1);
        nodes.push_back(n);
    }
    // setup the links
    for(int i = 0; i < node_generator->node_infos.size(); ++i) {
        // get them both in here
        auto node_info = node_generator->node_infos[i];
        auto node = nodes[i];
        //
        assert(node_info->links.size() > 1);
        int position = 0;
        for(int j = 1;  j < CHANNELS; ++j) {
            auto link = node_info->links[position];
            int channel = link.channel;
            handle_t partner_handle = link.handle;
            // actual partner
            if(channel == j) {
                long array_position = long(partner_handle) - 1;
                auto link_p = nodes[array_position]->populations.begin() + channel;
                node->neighbors.push_back(link_p);
                ++position;
            }
            // set to self
            else {
                channel = j;
                auto link_p = nodes[i]->populations.begin() + channel;
                node->neighbors.push_back(link_p);
            }
        }
    }
    // set up the markers
    handle_t maker_handles = 0;
    for (auto m : node_generator->markers->marker_points) {
        auto n = new marker(++maker_handles,*m);
        markers.push_back(n);
    }
    // set up marker pkh
    for(auto m : markers) {
        markers_pkh.fill_key(m->handle,m->position);
    }
    // set up sim parameters
    parameters.mean_marker_distance = original_markers->return_marker_distance();

}

void ibmSimulation::run(int current_step) {
    // offset control
    offset_sim = ((current_step +1) & 0x1) * 9;
    offset_node = (current_step & 0x1) * 9;
    for(auto n : nodes) {
        // todo where does the force actually come into the picture?!
        // shorthands
        auto force = n->forces;
        vector_t f = {0,0};
        std::vector<handle_t> markers_in_vicinity;
        auto population = n->populations;
        if(n->boundary_type == IBM) {
            // find the markers in the vicinity
            markers_in_vicinity = markers_pkh.ranging_key_translation(n->position,parameters.ibm_range);
            // ibm contribution to the force
            f = aggregate_force(&markers_in_vicinity,&n->position);
        }
        // calculate the macro values and do all the lbm stuff
        // macro
        auto [rho, ux, uy]  = calculate_macro_force(&population,&force);
        n->velocity = {ux, uy};
        // collision
        collision(&population,rho,ux,uy);
        // frogging
        if(n->boundary_type == NO_BOUNDARY) {
            f += calculate_rotation_force(&n->position,&n->velocity);
        }
        forcing_term(n,&f);
        // streaming
        streaming(&population,&n->neighbors);
        // distribute the new velocity to the markers
        if(n->boundary_type == IBM) {
            distribute_velocity(&markers_in_vicinity,&n->position,&n->velocity);
        }
    }
    // we propagate the markers forward and calculate the new forces based on the position of the marker
    propagate_calculate_force_marker();
}

void ibmSimulation::get_data(bool write_to_file) {
    /// flowfields
    flowfield_t ux;
    flowfield_t uy;
    flowfield_t rho;
    long size_x = size.x;
    long size_y = size.y;
    ux.resize(size_x,size_y);
    uy.resize(size_x,size_y);
    rho.resize(size_x,size_y);
    // make sure to get the correct one
    for(int i = 0; i < nodes.size(); ++i) {
        auto node = nodes[i];
        auto [rho_local,ux_local, uy_local] = calculate_macro_force(&node->populations,&node->forces);
        ux(int(node->position(0)),int(node->position(1))) = ux_local;
        uy(int(node->position(0)),int(node->position(1))) = uy_local;
        rho(int(node->position(0)),int(node->position(1))) = rho_local;
    }
    // write to a file otherwise useless
    write_flowfield_data(&ux, "ux_data_file",write_to_file);
    write_flowfield_data(&uy, "uy_data_file",write_to_file);
    write_flowfield_data(&rho, "rho_data_file",write_to_file);
}

void ibmSimulation::delete_containers() {
    for(auto n : nodes) {
        delete n;
    }
    for(auto m : markers) {
        delete m;
    }
    nodes.clear();
    markers.clear();
}

void ibmSimulation::test_propagate_markers() {
    propagate_calculate_force_marker();
}

double ibmSimulation::test_kernel_function(point_t *p) {
    return (*kernel_function)(p);
}

double ibmSimulation::test_kernel_function_call(point_t* p) {
    return kernel_function_call(p);
}

std::tuple<double,double,double> ibmSimulation::test_macro(array_t *a) {
    return calculate_macro(a);
}

void ibmSimulation::test_propagate_velocity() {
    for(auto n : nodes) {
        std::vector<handle_t> markers_in_vicinity;
        if(n->boundary_type == IBM) {
            // find the markers in the vicinity
            markers_in_vicinity = markers_pkh.ranging_key_translation(n->position,parameters.ibm_range);
            // ibm contribution to the velocity
            distribute_velocity(&markers_in_vicinity,&n->position,&n->velocity);
        }
    }
}