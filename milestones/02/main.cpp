/**
 * todos:
 * canvas scheme for initialization
 * grouping based on Ns
 * save initialization somewhere just do it once!!
 * incorporate more boundary methods
 * parallelize for equal nodes in a true NL structure
 * add doxygen docu, recheck old projects on a common standard
 */
#include "simulation.h"
#include <chrono>
#include <iostream>

// standard sliding lid
int main(int argc, char *argv[]) {
    auto start = std::chrono::high_resolution_clock::now();
    int size = 302;
    int steps = 10000;
    point_t p = {size,size};
    boundaryPointConstructor boundaries(p);
    boundaries.init_sliding_lid();
    // init the sim runner
    nodeGenerator gen(&boundaries);
    //gen.set_no_ordering();
    gen.init();
    // gen.set_no_ordering();
    simulation sim(&boundaries,&gen);
    sim.init();
    // init sim parameters
    double re = 1000;
    double base_length = size - 2;
    simulation_parameters params;
    params.u_wall = 0.1;
    params.relaxation = (2*re)/(6*base_length*params.u_wall+re);
    sim.set_simulation_parameters(params);
    // run sim
    for(int i = 0; i < steps; ++i) {
        sim.run();
        if(i % 1000 == 0) {
            std::cout << "Step: " << i << std::endl;
        }
    }
    sim.get_data(true,p);
    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds >(stop - start);
    std::cout << "Took " <<duration.count()<< " ms" << std::endl;
    return 0;
}
