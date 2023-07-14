#include "simulation.h"
#include "image_converter.h"
#include <iostream>
#include <chrono>
/**
 * 10 main, imb integration, we simulate a quader of a fluid rotating, boarders are off.
 * @param argc
 * @param argv
 * @return
 */
int main(int argc, char *argv[]) {
    point_t starter = {13.5,13.5};
    straight_t input;
    straightGenerator sg;
    int steps = 30;
    long canvas_size = 100;
    double side_length = 69; // with a distance of 0.75 we should get 80 markers
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
    ng.init_surface(canvas_size,2);
    ng.visualize_2D_nodes();
    std::cout << ng.node_infos.size() << std::endl;
    simulation_parameters params;
    params.relaxation = 0.8;
    point_t dk = {0,0};
    // max rotation is 7.5e-3
    vector_t sizes = {canvas_size,canvas_size};
    goaForce rot(dk,sizes,1e-5);
    forcedSimulation sim(&ng, &rot,ng.markers,sizes);;
    sim.init();
    for(int i = 0; i < steps; ++i) {
        if (i % 1000 == 0) {
            std::cout << "Step: " << i << std::endl;
        }
        sim.run_ibm(i);
    }
    sim.get_data(true);
    if(ng.straight_surfaces != nullptr)
        ng.straight_surfaces->write_out_surface();
    return 0;
}