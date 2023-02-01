#include "simulation.h"
#include <iostream>
#include <chrono>

int main(int argc, char *argv[]) {
    auto start = std::chrono::high_resolution_clock::now();
    int steps = 1000;
    unsigned int size = 32;
    unsigned int sub_size = 22;
    point_t c = {size,size};
    point_t p = {sub_size,sub_size+2};
    boundaryPointConstructor boundaries(p);
    // boundaries.init_sliding_lid_side_chopped({20,10},30);
    boundaries.init_sliding_lid_inner({3,2},{10,10},{10,13});
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
    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::seconds>(stop - start);
    std::cout << "Took " <<duration.count()<< "s" << std::endl;
    return 0;
}