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
#include <ctime>
#include <iostream>

int main(int argc, char *argv[]) {
    int size = 102;
    int steps = 10000;
    point_t p = {size,size};
    boundaryPointConstructor boundaries(p);
    boundaries.init_sliding_lid();
    // init the sim runner
    simulation sim(&boundaries);
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
    return 0;
}
