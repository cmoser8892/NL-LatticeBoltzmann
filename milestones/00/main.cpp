
/*
 * vallgrind call
valgrind --tool=callgrind --dump-instr=yes (p)
 */

#include <iostream>
#include "helper_functions.h"
#include "functions.h"

/**
 * Test main.
 * @return
 */
int main(int argc, char *argv[]) {
    double range = 0;
    double delta_range = 1;
    if(argc > 1) {
        range = std::stod(argv[1]);
        delta_range = std::stod(argv[2]);
    }
    std::cout << kernel_3(range,delta_range);
    return 0;
}