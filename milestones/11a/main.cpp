#include <iostream>
#include <chrono>
#include "lbm_simulation.h"
#include "simulation.h"
#include "one_step_simulation.h"

int main(int argc, char *argv[]) {
    point_t starter = {10,10};
    int steps = 3000;
    long canvas_size = 100;
    double side_length = 60; // with a distance of 0.75 we should get 80 markers
    point_t size = {canvas_size, canvas_size};
    point_t sinner = {side_length,side_length};
    boundaryPointConstructor boundaries(size);
    boundaries.init_quader(starter, sinner);
    nodeGenerator ng(&boundaries);
    ng.init_fused(canvas_size);
    ng.visualize_2D_nodes();
    simulation_parameters params;
    params.relaxation = 0.5;
    params.ibm_range = 2;
    point_t dk = {0,0};
    // max rotation is 7.5e-3
    goaForce rot(dk,size,1e-3);
    optimizedSimulation sim(&boundaries,&ng, &rot);
    sim.init();
    for(int i = 0; i < steps; ++i) {
        if(i % 1000 == 0) {
            std::cout << "Step: " << i << std::endl;
        }
        sim.forcing_run(i);
    }
    sim.get_data(true);
    if(ng.straight_surfaces != nullptr)
        ng.straight_surfaces->write_out_surface();
}
