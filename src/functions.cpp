#include "functions.h"
#include <cmath>
/// globals
double relaxation = 0.5;
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
// create a copy array copy values into data (asynchron) copy from old array to new one
void streaming_step1(node* node) {
    if(node->neighbors.empty()) {
        return;
    }
    for( int i = 1; i < node->data.size(); ++i) {
        if(node->neighbors.at(i-1) != nullptr) {
            node->neighbors.at(i-1)->copy(i) = node->data(i);
        }
    }
}

void streaming_step2(node* node) {
    node->data = node->copy;
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

void moving_wall(node * node,int side_pos,double uw) {
    /// only for y + i hate it couse very manual
    // we look for boundary nodes and then see where they want to go maximum would be 3
    // maybe lock equation to the node?! more function pointers todo
    if (node->node_type == BOUNDARY) {
        if(node->position.y() == side_pos) {
            int count = 0;
            for( int i = 1; i < node->data.size(); ++i) {
                if(node->neighbors.at(i-1) != nullptr) {
                    if(count == 0) {
                        ++count;
                    }
                    else if(count == 1) {
                        ++count;
                        node->neighbors.at(i-1)->copy(i) -= 1.0/6*uw;
                    }
                    else if(count == 2) {
                        node->neighbors.at(i-1)->copy(i) += 1.0/6*uw;
                    }
                }
            }
        }
    }
}

void write_ux(node* node, flowfield_t* ux) {
    ux->operator()(int(node->position(0)),int(node->position(1))) = node->u(0);
}

void write_uy(node* node, flowfield_t * uy) {
    uy->operator()(int(node->position(0)),int(node->position(1))) = node->u(1);
}