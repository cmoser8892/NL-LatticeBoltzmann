#include "one_step_simulation.h"
#include <iostream>
#include <chrono>

/**
 * 06b main,special pretest version that doesnt do bb moving to investigate the speedups of direct array access.
 * @param argc
 * @param argv
 * @return
 */
int main(int argc, char *argv[]) {
    auto start = std::chrono::high_resolution_clock::now();
    // init variant
    int steps = 5;
    unsigned int size = 200;
    point_t c = {size,size};
    boundaryPointConstructor boundaries(c);
    // boundaries.init_sliding_lid_side_chopped({20,10},30);
    boundaries.init_sliding_lid();
    nodeGenerator gen(&boundaries);
    gen.init_fused(size);
    // init sim parameters
    double re = 1000;
    double base_length = size - 2;
    simulation_parameters params;
    params.u_wall = 0.1;
    params.relaxation = (2*re)/(6*base_length*params.u_wall+re);
    optimizedSimulation sim(&boundaries,&gen);
    sim.set_simulation_parameters(params);
    sim.init();
    sim.init_sub_array();
    // run sim
    for(int i = 0; i < steps; ++i) {
        if(i % 1000 == 0) {
            std::cout << "Step: " << i << std::endl;
        }
        sim.run_sub_array(i);
    }
    sim.get_data(true,c);
    if(gen.straight_surfaces != nullptr)
        gen.straight_surfaces->write_out_surface();
    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);
    std::cout << "Took " <<duration.count()<< "ms" << std::endl;
    return 0;
}