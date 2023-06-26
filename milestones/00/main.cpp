
/*
 * vallgrind call
valgrind --tool=callgrind --dump-instr=yes (p)
 */

#include <iostream>
#include "helper_functions.h"

/**
 * Test main.
 * @return
 */
int main() {
    double size = 300;
    double re = 1000;
    double base_length = size - 2;
    simulation_parameters params;
    params.u_wall = 0.1;
    params.relaxation = (2*re)/(6*base_length*params.u_wall+re);
    std::cout << params.relaxation << std::endl;
    // recalculate the actual reynolds number of the slinding lid with a hole
    double viscocity = 1.0/3 * (1/params.relaxation -0.5);
    double sub_length = 220;
    re = (sub_length * params.u_wall)/viscocity;
    std::cout << re << std::endl;
}