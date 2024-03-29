#include <iostream>
#include <chrono>
#include "lbm_simulation.h"

/**
 * 11b main, the populations in the nodes on the edge of the simulation domain have a bounce back enable.
 * @result No real difference to 11, when working with a circular force, disappointing.
 * @param argc
 * @param argv
 * @return
 */
int main(int argc, char *argv[]) {
    point_t starter = {13,13};
    straight_t input;
    straightGenerator sg;
    int steps = 10000;
    long canvas_size = 100;
    double marker_distance = 0.5;
    double side_length = 60; // with a distance of 0.75 we should get 80 markers
    kernelType_t kernel = KERNEL_C;
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
    /// todo main diff between 11 and 11b
    ng.init_surface_return(canvas_size,kernel,marker_distance);
    ng.visualize_2D_nodes();
    std::cout << ng.node_infos.size() << std::endl;
    simulation_parameters params;
    params.relaxation = 0.5;
    params.ibm_range = kernel_id_to_lattice_search(kernel);
    params.kernel_in_use = kernel;
    params.k = 1;
    point_t dk = {0,0};
    // max rotation is 7.5e-3
    vector_t sizes = {canvas_size,canvas_size};
    goaForce rot(dk,sizes,1e-3);
    ibmSimulation sim(&ng, &rot,ng.markers,sizes);
    sim.set_simulation_parameters(params);
    sim.init();
    // put the markers in the watchdog
    rhoWatchdog dog(0.1,sizes);
    markerWatchdog marker_watch(0.1);
    long marker_check = 0;
    for(auto m : sim.markers) {
        marker_watch.init_marker(m->position);
    }
    // run the sim
    for(int i = 0; i < steps; ++i) {
        if(i % 1000 == 0) {
            std::cout << "Step: " << i << std::endl;
        }
        sim.run(i);
        for(auto m : sim.markers) {
            if(marker_watch.check(m->position,m->handle) == true) {
                marker_check++;
                if(1) {
                    std::cout << m->position.x() << ", " << m->position.y() << "//"
                              << m->original_position.x() << ", " << m->original_position.y() << std::endl;
                }
            }
        }
        // watchdog part
        if(1) {
            for(auto n : sim.nodes) {
                if(dog.check_force(n,i)) {
                }
            }
        }
    }

    sim.get_data(true);
    if(ng.straight_surfaces != nullptr)
        ng.straight_surfaces->write_out_surface();
    std::cout << std::endl;
    std::cout << std::endl;
    std::cout << std::endl;
    for(auto m : sim.markers) {
        std::cout << m->position.x() << ", " << m->position.y() << "//"
                  << m->original_position.x() << ", " << m->original_position.y() << std::endl;
    }
    std::cout << marker_check << std::endl;
}
