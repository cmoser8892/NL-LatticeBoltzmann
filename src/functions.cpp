#include "functions.h"
#include <cmath>
/// globals
double relaxation = 0.5;
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
    return_array(0) = 2.0/9 * node->rho * (2-3*uu);
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

void streaming(node* node) {
    // prob the hard method
}

void collision(node* node) {
    node->data -= relaxation*(node->data - equilibrium(node));
}

void macro( node* node) {
    node->rho = node->data.sum();
    node->u(0) = ((node->data(1)+node->data(5)+node->data(8))-
                  (node->data(3)+node->data(6) +node->data(7)))/node->rho;
    node->u(1) = ((node->data(2)+node->data(5)+node->data(6))-
                  (node->data(4)+node->data(7) +node->data(8)))/node->rho;
}