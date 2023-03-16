#include "functions.h"
#include <cmath>
#include <iostream>

/// globals
matrix_t velocity_set = {{0,1,0,-1,0 ,1,-1,-1, 1},
                         {0,0,1,0 ,-1,1, 1,-1,-1}};

matrix_t weights = {{4.0/9,1.0/9,1.0/9,1.0/9,1.0/9,1.0/36,1.0/36,1.0/36,1.0/36}};

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
 * @fn void one_step_macro_collision(oNode* node, double relaxation)
 * @brief one step for all the calculations necessary
 * @param node
 * @param relaxation
 */
void one_step_macro_collision(oNode* node, double relaxation) {
    int o = node->offset;
    // macro calc
    auto p = node->populations.begin() + o;
    // macro part
    // rho
    const double rho = (p + 0).operator*() +
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
    // collision
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
 * @fn double calculate_ux(oNode* node)
 * @brief calculate ux
 * @param node
 * @return
 */
double calculate_ux(oNode* node) {
    int o = node->offset;
    auto p = node->populations.begin() + o;
    double rho = calculate_rho(node);
    double ux = (((p + 1).operator*() +
                  (p + 5).operator*() +
                  (p + 8).operator*())-
                 ((p + 3).operator*() +
                  (p + 6).operator*()+
                  (p + 7).operator*()));
    return ux/rho;
}

/**
 * @fn double calculate_uy(oNode* node)
 * @brief calculates uy
 * @param node
 * @return
 */
double calculate_uy(oNode* node) {
    int o = node->offset;
    auto p = node->populations.begin() + o;
    double rho = calculate_rho(node);
    double uy = (((p + 2).operator*() +
                  (p + 5).operator*() +
                  (p + 6).operator*())-
                 ((p + 4).operator*() +
                  (p + 7).operator*()+
                  (p + 8).operator*()));
    return uy/rho;
}

/**
 * @fn double calculate_rho(oNode* node)
 * @brief calculates rho
 * @param node
 * @return
 */
double calculate_rho(oNode* node) {
    int o = node->offset;
    auto p = node->populations.begin() + o;
    double rho = (p + 0).operator*() +
                 (p + 1).operator*() +
                 (p + 2).operator*() +
                 (p + 3).operator*() +
                 (p + 4).operator*() +
                 (p + 5).operator*() +
                 (p + 6).operator*() +
                 (p + 7).operator*() +
                 (p + 8).operator*();
    return rho;
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

/**
 * @fn rhoWatchdog::rhoWatchdog(double s,point_t size
 * @brief constructor for the rho_watchdog
 * @param s
 * @param size
 */
rhoWatchdog::rhoWatchdog(double s,point_t size) :sensitivity(s) {
    rho.setOnes(long(size.x()),long(size.y()));
}

/**
 * @fn bool rhoWatchdog::check(node *n,int step)
 * @brief performs a watchdog check of the history of the rho value, aka compares it to the previous one
 * @param n
 * @param step
 * @return
 */
bool rhoWatchdog::check(node *n,int step) {
    double rho_old = rho(int(n->position(0)),int(n->position(1)));
    bool return_value = false;
    if((abs(rho_old-n->rho)) >= (abs(rho_old*sensitivity))) {
        std::cerr << "Rho-diviation at " << step << std::endl;
        std::cerr << "Position: " << n->position.x()
                  << " ," << n->position.y()
                  << std::endl;
        std::cerr << "Rho previous: " << rho_old << std::endl;
        std::cerr << "Rho now: " << n->rho << std::endl;
        std::cerr << std::endl;
        return_value = true;
    }
    rho(int(n->position(0)),int(n->position(1))) = n->rho;
    return return_value;
}