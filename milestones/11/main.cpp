#include <iostream>
#include <chrono>
#include "lbm_simulation.h"

int main(int argc, char *argv[]) {
    point_t starter = {13,13};
    straight_t input;
    straightGenerator sg;
    int steps = 300;
    long canvas_size = 100;
    double side_length = 61; // with a distance of 0.75 we should get 80 markers
    double ibm_distance = 2;
    // we put in a quader
    input.point = starter;
    input.direction = {0,1};
    input.max_t = side_length;
    sg.add_surface(input);
    input.point += input.direction * side_length;
    input.direction = {1,0};
    input.max_t = side_length;
    sg.add_surface(input);
    input.point += input.direction * side_length;
    input.direction = {0,-1};
    input.max_t = side_length;
    sg.add_surface(input);
    input.point += input.direction * side_length;
    input.direction = {-1,0};
    input.max_t = side_length;
    sg.add_surface(input);
    sg.surface_mass_center();
    nodeGenerator ng(&sg);
    ng.init_surface(canvas_size,ibm_distance);
    ng.visualize_2D_nodes();
    simulation_parameters params;
    params.relaxation = 0.5;
    params.ibm_range = ibm_distance;
    params.k = 1;
    point_t dk = {0,0};
    // max rotation is 7.5e-3
    vector_t sizes = {canvas_size,canvas_size};
    goaForce rot(dk,sizes,1e-3);
    ibmSimulation sim(&ng, &rot,ng.markers,sizes);
    sim.set_simulation_parameters(params);
    sim.init();
    for(int i = 0; i < steps; ++i) {
        if(i % 1000 == 0) {
            std::cout << "Step: " << i << std::endl;
        }
        sim.run(i);
    }
    sim.get_data(true);
    if(ng.straight_surfaces != nullptr)
        ng.straight_surfaces->write_out_surface();
}
