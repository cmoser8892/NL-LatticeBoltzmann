#include "functions.h"
#include <cmath>
#include <iostream>

/// globals
matrix_t velocity_set = {{0,1,0,-1,0 ,1,-1,-1, 1},
                         {0,0,1,0 ,-1,1, 1,-1,-1}};

/// implementation is for d2q9
array_t equilibrium(node* node) {
    array_t return_array;
    return_array.setZero(node->data.size());
    array_t pow = node->u.pow(2);
    double uu = 3*(pow.sum());
    double six_ux = 6 * node->u(0);
    double six_uy = 6 * node->u(1);
    double nine_uxx = 9 * pow(0);
    double nine_uyy = 9 * pow(1);
    double three_ux_uy_p = 3 * ( node->u(0) + node->u(1));
    double three_ux_uy_m = 3 * ( node->u(0) - node->u(1));
    double nine_uxuy = 9 * ( node->u(0) * node->u(1));
    return_array(0) = 2.0/9 * node->rho * (2- uu);
    return_array(1) = 1.0/18 * node->rho * (2 + six_ux + nine_uxx - uu);
    return_array(2) = 1.0/18 * node->rho * (2 + six_uy + nine_uyy - uu);
    return_array(3) = 1.0/18 * node->rho * (2 - six_ux + nine_uyy - uu);
    return_array(4) = 1.0/18 * node->rho * (2 - six_uy + nine_uyy - uu);
    return_array(5) = 1.0/36 * node->rho * (1 + three_ux_uy_p + nine_uxuy + uu);
    return_array(6) = 1.0/36 * node->rho * (1 - three_ux_uy_m - nine_uxuy + uu);
    return_array(7) = 1.0/36 * node->rho * (1 - three_ux_uy_p + nine_uxuy + uu);
    return_array(8) = 1.0/36 * node->rho * (1 + three_ux_uy_m - nine_uxuy + uu);
    return return_array;
}

void collision(node* node, double relaxation) {
    node->data -= relaxation * (node->data - equilibrium(node));
}

void macro( node* node) {
    node->rho = node->data.sum(); // problem if 0!!
    node->u(0) = ((node->data(1)+node->data(5)+node->data(8))-
                  (node->data(3)+node->data(6) +node->data(7)));
    node->u(0) = node->u(0)/node->rho;
    node->u(1) = ((node->data(2)+node->data(5)+node->data(6))-
                  (node->data(4)+node->data(7) +node->data(8)))
                 /node->rho;
    node->u(1) = node->u(1)/node->rho;
}

// write the ux component of the flowfield
void write_ux(node* node, flowfield_t* ux) {
    // dont ask this looks ugly
    ux->operator()(int(node->position(0)),int(node->position(1))) = node->u(0);
}

// writes the uy_component of a flowfield
void write_uy(node* node, flowfield_t * uy) {
    // dont ask this looks ugly
    uy->operator()(int(node->position(0)),int(node->position(1))) = node->u(1);
}

void write_rho(node* node, flowfield_t * rho) {
    // dont ask this looks ugly
    rho->operator()(int(node->position(0)),int(node->position(1))) = node->rho;
}

void debug_node(node* node, bool printing) {
    // print out the data values and calculate the density
    if(printing) {
        std::cout << "Data" << std::endl;
        std::cout << node->data << std::endl;
        std::cout << "rho" << std::endl;
        std::cout << node->rho << std::endl<< std::endl;
    }
}

double bb_switch_channel(int from_channel, double uw) {
    // incomplete all around implementation only for top side aka channels 7 and 8
    double return_value = 0;
    switch(from_channel) {
    case 5:
        return_value = -1.0/6 * uw;
        break;
    case 6:
        return_value = 1.0/6 * uw;
        break;
    default:
        // nop
        break;
    }
    return return_value;
}