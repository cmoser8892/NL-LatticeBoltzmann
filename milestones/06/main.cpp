#include "simulation.h"
#include <iostream>
#include <chrono>

/// fused variant of m3
int main(int argc, char *argv[]) {
    auto start = std::chrono::high_resolution_clock::now();
    int steps = 10000;
    unsigned int size = 302;
    unsigned int sub_size = 202;
    point_t c = {size,size};
    point_t p = {sub_size,sub_size+20};
    boundaryPointConstructor boundaries(p);
    // boundaries.init_sliding_lid_side_chopped({20,10},30);
    boundaries.init_sliding_lid_inner({10,20},{34,45},{49,52});
    nodeGenerator gen(&boundaries);
    gen.init_fused(size);
    // init sim parameters
    double re = 1000;
    double base_length = size - 2;
    simulation_parameters params;
    params.u_wall = 0.1;
    params.relaxation = (2*re)/(6*base_length*params.u_wall+re);
    oSimu sim(&boundaries,&gen);
    sim.set_simulation_parameters(params);
    sim.init();
    // run sim
    for(int i = 0; i < steps; ++i) {
        sim.run(i);
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