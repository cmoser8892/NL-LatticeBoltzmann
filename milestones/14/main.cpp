
/*
 * vallgrind call
valgrind --tool=callgrind --dump-instr=yes (p)
 */

#include "lbm_simulation.h"
#include "helper_functions.h"
#include "drawn_image_surface.h"
#include "nodeGenerator.h"

/**
 * 14 pflow with a force acting on the nodes with a costume surface
 * @return
 */
int main() {
    // node generator variables
    long canvas_size = 150;
    vector_t draw_size = {canvas_size-1,canvas_size-1};
    double marker_distance = 0.5;
    bool file_write = true;
    kernelType_t kernel = KERNEL_C;
    // introduce a surface
    straight_t input;
    straightGenerator lines;
    // set the input and introduce
    input.point = {0,27.1};
    input.direction = {1,0};
    input.max_t = 149;
    input.type = IBM;
    lines.add_surface(input);
    input.point = {149,82.8};
    input.direction = {-1,0};
    input.max_t = 149;
    input.type = IBM;
    lines.add_surface(input);
    // introduce the containing lines will have to be manually found !
    input.point = {0,27.1};
    input.direction = {0,1};
    input.max_t = 82.8-27.1;
    input.type = PERIODIC;
    lines.add_surface(input);
    input.point = {149,27.1};
    input.direction = {0,1};
    input.max_t = 82.8-27.1;
    input.type = PERIODIC;
    lines.add_surface(input);
    // setup the node generator
    nodeGenerator ng(&lines);
    ng.init_surface_return(canvas_size,kernel,marker_distance);
    // write out all the boundary types found
    ng.write_out_nodes(IBM_INNER, file_write);
    ng.write_out_nodes(IBM_OUTER, file_write);
    ng.write_out_nodes(NO_BOUNDARY, file_write);
    // setup params for sim
    int steps = 15000;
    simulation_parameters params;
    params.relaxation = 0.5;
    params.ibm_range = kernel_id_to_lattice_search(kernel);
    params.kernel_in_use = kernel;
    params.k = 1;
    vector_t sizes = {canvas_size,canvas_size};
    ibmSimulation sim(&ng, nullptr,ng.markers,sizes);
    sim.set_simulation_parameters(params);
    sim.init();
    // run
    for(int i = 0; i <steps; ++i) {
        if(i % 100 == 0) {
            std::cout << "Step: " << i << std::endl;
        }
        sim.run_snake(i);
    }
    //
    sim.get_data(true);
    if(ng.straight_surfaces != nullptr)
        ng.straight_surfaces->write_out_surface();
    // end
    return 0;
}
