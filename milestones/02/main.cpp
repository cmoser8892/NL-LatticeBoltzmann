#include "two_step_simulation.h"
#include <chrono>
#include <iostream>

/**
 * Demonstration of a basic sliding lid.
 * @param argc
 * @param argv
 * @return
 */
int main(int argc, char *argv[]) {
    auto start = std::chrono::high_resolution_clock::now();
    int size = 302;
    int steps = 10000;
    point_t p = {size,size};
    boundaryPointConstructor boundaries(p);
    boundaries.init_sliding_lid();
    // init the sim runner
    nodeGenerator gen(&boundaries);
    gen.init();
    // init sim parameters
    double re = 1000;
    double base_length = size - 2;
    simulation_parameters params;
    params.u_wall = 0.1;
    params.relaxation = (2*re)/(6*base_length*params.u_wall+re);
    basicSimulation sim(&boundaries,&gen);
    sim.set_simulation_parameters(params);
    sim.init();
    // run sim
    for(int i = 0; i < steps; ++i) {
        bool once = false;
        if(i % 1000 == 0) {
            std::cout << "Step: " << i << std::endl;
            // sim.get_data(false,p);
        }
        sim.run();
        /*
        // watchdog part
        for(auto n : sim.nodes) {
            if(dog->check(n,i)) {
                once = true;
            }
        }
        if(once) {
            sim.get_data(false,p);
        }
         */
    }
    sim.get_data(true,p);
    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds >(stop - start);
    std::cout << "Took " <<duration.count()<< " ms" << std::endl;
    return 0;
}
