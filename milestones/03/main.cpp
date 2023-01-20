#include "simulation.h"
#include <iostream>

int main(int argc, char *argv[]) {
    int steps = 10000;
    unsigned int size = 132;
    unsigned int sub_size = 102;
    point_t c = {size,size};
    point_t p = {sub_size,sub_size};
    boundaryPointConstructor boundaries(p);
    boundaries.init_chopped_sliding_lid({20,10},0);
    nodeGenerator gen(&boundaries);
    gen.init(size);
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
    sim.get_data(true,c);

    return 0;
}
