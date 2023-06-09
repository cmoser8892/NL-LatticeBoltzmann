#include "functions.h"
#include <iostream>

/// implementation is for d2q9
/**
 * @fn array_t equilibrium(node* node)
 * @brief implementation of the equilibrium function
 * @param node
 * @return
 */
array_t equilibrium(node* node) {
    // this optimizes badly
    array_t return_array;
    return_array.setZero(CHANNELS);
    // convenience programming
    double ux = node->u(0);
    double uy = node->u(1);
    double rho = node->rho;
    // unroll the collision function
    return_array(0) = weights.col(0).x()*rho*(1- 1.5*(ux*ux +uy*uy));
    return_array(1) = weights.col(1).x()*rho*(1+ 3*ux+ 4.5*ux*ux- 1.5*(ux*ux +uy*uy));
    return_array(2) = weights.col(2).x()*rho*(1+ 3*uy+ 4.5*uy*uy- 1.5*(ux*ux +uy*uy));
    return_array(3) = weights.col(3).x()*rho*(1- 3*ux+ 4.5*ux*ux- 1.5*(ux*ux +uy*uy));
    return_array(4) = weights.col(4).x()*rho*(1- 3*uy+ 4.5*uy*uy- 1.5*(ux*ux +uy*uy));
    return_array(5) = weights.col(5).x()*rho*(1+ 3*ux+ 3*uy+ 9*ux*uy+ 3*(ux*ux +uy*uy));
    return_array(6) = weights.col(6).x()*rho*(1- 3*ux+ 3*uy- 9*ux*uy+ 3*(ux*ux +uy*uy));
    return_array(7) = weights.col(7).x()*rho*(1- 3*ux- 3*uy+ 9*ux*uy+ 3*(ux*ux +uy*uy));
    return_array(8) = weights.col(8).x()*rho*(1+ 3*ux- 3*uy- 9*ux*uy+ 3*(ux*ux +uy*uy));
    return return_array;
}

/**
 * @fn array_t equilibrium_general(node* node)
 * @brief original version of the equiliburm function as implemented in lookup
 * @param node
 * @return
 */
array_t equilibrium_general(node* node) {
    array_t return_array;
    return_array.setZero(CHANNELS);
    double ux = node->u(0);
    double uy = node->u(1);
    double rho = node->rho;
    for (int i = 0; i < CHANNELS; ++i) {
        double cx = velocity_set.col(i).x();
        double cy = velocity_set.col(i).y();
        double w = weights.col(i).x();
        return_array(i) = w * rho * calculate_later_equilibrium(cx, cy, ux, uy);
    }
    return return_array;
}

/**
 * @fn double calculate_later_equilibrium(double cx, double cy, double ux, double uy)
 * @brief calculates the later part of the equilbirum function
 * @param cx
 * @param cy
 * @param ux
 * @param uy
 * @return
 */
double calculate_later_equilibrium(double cx, double cy, double ux, double uy) {
    return (1 + 3*cx*ux + 3*cy*uy + 4.5*cx*cx*ux*ux + 9*cx*ux*cy*uy + 4.5*cy*cy*uy*uy - 1.5*ux*ux - 1.5*uy*uy);
}

/**
 * @fn array_t equilibrium_2d(double ux, double uy, double rho)
 * @brief calculates the equilibrium function based on macro values
 * @param ux
 * @param uy
 * @param rho
 * @return
 */
array_t equilibrium_2d(double ux, double uy, double rho) {
    array_t return_array;
    return_array.setZero(CHANNELS);
    return_array(0) = weights.col(0).x()*rho*(1- 1.5*(ux*ux +uy*uy));
    return_array(1) = weights.col(1).x()*rho*(1+ 3*ux+ 4.5*ux*ux- 1.5*(ux*ux +uy*uy));
    return_array(2) = weights.col(2).x()*rho*(1+ 3*uy+ 4.5*uy*uy- 1.5*(ux*ux +uy*uy));
    return_array(3) = weights.col(3).x()*rho*(1- 3*ux+ 4.5*ux*ux- 1.5*(ux*ux +uy*uy));
    return_array(4) = weights.col(4).x()*rho*(1- 3*uy+ 4.5*uy*uy- 1.5*(ux*ux +uy*uy));
    return_array(5) = weights.col(5).x()*rho*(1+ 3*ux+ 3*uy+ 9*ux*uy+ 3*(ux*ux +uy*uy));
    return_array(6) = weights.col(6).x()*rho*(1- 3*ux+ 3*uy- 9*ux*uy+ 3*(ux*ux +uy*uy));
    return_array(7) = weights.col(7).x()*rho*(1- 3*ux- 3*uy+ 9*ux*uy+ 3*(ux*ux +uy*uy));
    return_array(8) = weights.col(8).x()*rho*(1+ 3*ux- 3*uy- 9*ux*uy+ 3*(ux*ux +uy*uy));
    return return_array;
}
/**
 * @fn void collision(node* node, double relaxation)
 * @brief implementation of the collision function
 * @param node
 * @param relaxation
 */
void collision(node* node, double relaxation) {
    node->population_even -= relaxation * (node->population_even - equilibrium_general(node));
}

/**
 * @fn void macro( node* node)
 * @brief calculates the macroscopic values for the equilibrium function
 * @param node
 */
void macro( node* node) {
    node->rho = node->population_even.sum(); // problem if 0!!
    node->u(0) = ((node->population_even(1)+node->population_even(5)+node->population_even(8))-
                  (node->population_even(3)+node->population_even(6) +node->population_even(7)));
    node->u(1) = ((node->population_even(2)+node->population_even(5)+node->population_even(6))-
                  (node->population_even(4)+node->population_even(7) +node->population_even(8)));
    node->u /= node->rho; // will produce nonsense if div through 0
}

/**
 * @fn void fused_macro( node* node)
 * @brief calculates the macroscopic values for the equilibrium function for the fused runs
 * @param node
 */
void fused_macro(node * node) {
    node->rho = node->current_population->sum();
    node->u(0) = ((node->current_population->operator()(1) +
                   node->current_population->operator()(5) +
                   node->current_population->operator()(8))-
                  (node->current_population->operator()(3) +
                   node->current_population->operator()(6) +
                   node->current_population->operator()(7)));

    node->u(1) = ((node->current_population->operator()(2) +
                   node->current_population->operator()(5) +
                   node->current_population->operator()(6))-
                  (node->current_population->operator()(4) +
                   node->current_population->operator()(7) +
                   node->current_population->operator()(8)));

    node->u /= node->rho;
}

/**
 * @fn void fused_collision(node* node, double relax)
 * @brief variant of the optimized collision
 * @param node
 * @param relax
 */
void fused_collision(node* node, double relax) {
    // convenience programming
    double ux = node->u(0);
    double uy = node->u(1);
    double rho = node->rho;
    // unroll the collision function also messy for a reason optimizes better :)
    node->current_population->operator()(0) -= relax * (node->current_population->operator()(0) - weights.col(0).x()*rho*(1- 1.5*(ux*ux +uy*uy)));
    node->current_population->operator()(1) -= relax * (node->current_population->operator()(1) - weights.col(1).x()*rho*(1+ 3*ux+ 4.5*ux*ux- 1.5*(ux*ux +uy*uy)));
    node->current_population->operator()(2) -= relax * (node->current_population->operator()(2) - weights.col(2).x()*rho*(1+ 3*uy+ 4.5*uy*uy- 1.5*(ux*ux +uy*uy)));
    node->current_population->operator()(3) -= relax * (node->current_population->operator()(3) - weights.col(3).x()*rho*(1- 3*ux+ 4.5*ux*ux- 1.5*(ux*ux +uy*uy)));
    node->current_population->operator()(4) -= relax * (node->current_population->operator()(4) - weights.col(4).x()*rho*(1- 3*uy+ 4.5*uy*uy- 1.5*(ux*ux +uy*uy)));
    node->current_population->operator()(5) -= relax * (node->current_population->operator()(5) - weights.col(5).x()*rho*(1+ 3*ux+ 3*uy+ 9*ux*uy+ 3*(ux*ux +uy*uy)));
    node->current_population->operator()(6) -= relax * (node->current_population->operator()(6) - weights.col(6).x()*rho*(1- 3*ux+ 3*uy- 9*ux*uy+ 3*(ux*ux +uy*uy)));
    node->current_population->operator()(7) -= relax * (node->current_population->operator()(7) - weights.col(7).x()*rho*(1- 3*ux- 3*uy+ 9*ux*uy+ 3*(ux*ux +uy*uy)));
    node->current_population->operator()(8) -= relax * (node->current_population->operator()(8) - weights.col(8).x()*rho*(1+ 3*ux- 3*uy- 9*ux*uy+ 3*(ux*ux +uy*uy)));
}

/**
 * @fn void write_ux(node* node, flowfield_t* ux)
 * @brief write the ux component of the flowfield
 * @param node
 * @param ux
 */
void write_ux(node* node, flowfield_t* ux) {
    // dont ask this looks ugly
    ux->operator()(int(node->position(0)),int(node->position(1))) = node->u(0);
}

/**
 * @fn void write_uy(node* node, flowfield_t * uy)
 * @brief writes the uy_component of a flowfield
 * @param node
 * @param uy
 */
void write_uy(node* node, flowfield_t * uy) {
    // dont ask this looks ugly
    uy->operator()(int(node->position(0)),int(node->position(1))) = node->u(1);
}

/**
 * @fn void write_rho(node* node, flowfield_t * rho)
 * @brief writes the rho component of a flow field
 * @param node
 * @param rho
 */
void write_rho(node* node, flowfield_t * rho) {
    // dont ask this looks ugly
    rho->operator()(int(node->position(0)),int(node->position(1))) = node->rho;
}

/**
 * @fn void debug_node(node* node, bool printing)
 * @brief print out stuff dont forget to set printing to true
 * @param node
 * @param printing
 */
void debug_node(node* node, bool printing) {
    // print out the population_even values and calculate the density
    if(printing) {
        std::cout << "Position" << std::endl;
        std::cout << node->position << std::endl;
        std::cout << "population_even" << std::endl;
        std::cout << node->population_even << std::endl;
        std::cout << "rho" << std::endl;
        std::cout << node->rho << std::endl<< std::endl;
    }
}

/**
 * @fn double bb_switch_channel(int from_channel, double uw)
 * @brief applies the moving part of the boundary to a channel
 * @param from_channel
 * @param uw
 * @return how much needs to be subtracted for moving boundaries
 */
double bb_switch_channel(int from_channel, double uw) {
    // incomplete all around implementation only for top side aka channels 7 and 8
    double return_value = 0;
    switch(from_channel) {
    case 5:
        return_value = -1.0/6 * uw;
        break;
    case 6:
        return_value = +1.0/6 * uw;
        break;
    default:
        // nop
        break;
    }
    return return_value;
}

/**
 * @fn int switch_link_dimensions(int link_channel)
 * @brief changes around the channel in a bounce back call
 * @param link_channel
 * @return return the right link channel
 */
int switch_link_dimensions(int link_channel) {
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

std::tuple<double, double ,double> calculate_force_macro_values(array_t p, array_t f) {
    double rho = 0;
    double ux = 0;
    double uy = 0;
    /// looped version to find errors in the other one
    //rho
    for(int i = 0; i < CHANNELS; ++i) {
        rho += (p[i] + 0.5*f[i]);
    }
    // velocity
    vector_t u;
    u.setZero();
    for(int i  = 0; i < CHANNELS; ++i) {
        vector_t v_set = velocity_set.col(i);
        u += p[i]*v_set;
    }
    // additional force term
    // add half of f_alpha, we have to recover that here prob better to just add it in post
    for(int i = 0; i < CHANNELS; ++i) {
        vector_t v_set = velocity_set.col(i);
        u += 1.0/2*f(i)*v_set;
    }
    return {rho,ux,uy};
}

std::tuple<double, double, double>  calculate_the_macro(array_t* a, int offset) {
    // return als a struct
    // calculate rho ux and uy
    int o = offset;
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

// todo implement me
void pressure_periodic_in(oNode* node, double rho_in) {
    int o = 0;
    // macro calc
    auto p = node->populations.begin() + o;
}

void pressure_periodic_out(oNode* node , double rho_out) {
    int o = 0;
    // macro calc
    auto p = node->populations.begin() + o;
    // array_t eq = equilibrium_2d(ux,uy,rho);
    // array_t eq_out = equilibrium_2d(ux,uy,rho_out);
    // pressure is not local :/
    for(int i = 0; i < CHANNELS; ++i) {
        (p+i).operator*() = 0;
    }
}
/*
// python stuff
def periodic_boundary_with_pressure_variations(grid,rho_in,rho_out):
    # get all the values
    rho, ux, uy = caluculate_real_values(grid)
    equilibrium = equilibrium_on_array(rho, ux, uy)
    ##########
    equilibrium_in = equilibrium_on_array(rho_in, ux[-2,:], uy[-2, :])
    # inlet 1,5,8
    grid[:, 0, :] = equilibrium_in + (grid[:, -2, :] - equilibrium[:, -2, :])

# outlet 3,6,7
equilibrium_out = equilibrium_on_array(rho_out, ux[1, :], uy[1, :])
# check for correct sizes
    grid[:, -1, :] = equilibrium_out + (grid[:, 1, :] - equilibrium[:, 1, :])
                                           */