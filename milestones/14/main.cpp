
/*
 * vallgrind call
valgrind --tool=callgrind --dump-instr=yes (p)
 */

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
    double ibm_distance = kernel_id_to_lattice_search(kernel);
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
    ng.init_surface_return(canvas_size,ibm_distance,marker_distance);
    // write out all the boundary types found
    ng.write_out_nodes(IBM_INNER, file_write);
    ng.write_out_nodes(IBM_OUTER, file_write);
    ng.write_out_nodes(NO_BOUNDARY, file_write);
    // write out the markers
    ng.markers->write_out_markers(file_write);
    std::cout << ng.markers->marker_points.size() << std::endl;
    lines.write_out_surface();
    // end
    return 0;
}
