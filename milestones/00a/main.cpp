
/*
 * vallgrind call
valgrind --tool=callgrind --dump-instr=yes (p)
 */

#include "helper_functions.h"
#include "drawn_image_surface.h"
#include "nodeGenerator.h"

/**
 * 00 main, testing stuff.
 * @note good rule of thumb is all 4 contours need to be selected
 * @return
 */
int main() {
    // node generator variables
    long canvas_size = 50;
    double marker_distance = 0.5;
    bool file_write = true;
    kernelType_t kernel = KERNEL_C;
    double ibm_distance = kernel_id_to_lattice_search(kernel);
    // Load the image
    straight_t input;
    straightGenerator sg;
    double side_length = 35; // with a distance of 0.75 we should get 80 markers
    point_t starter = {5,5};
    // we put in a quader
    input.point = starter;
    input.direction = {0,1};
    input.max_t = side_length;
    sg.add_surface(input);
    input.point += input.direction * side_length;
    vector_t v = {1,1};
    input.direction = v.normalized();
    input.max_t = 10;
    sg.add_surface(input);
    sg.surface_mass_center();
    // node genertor
    nodeGenerator ng(&sg);
    ng.init_surface_return(canvas_size,ibm_distance,marker_distance);
    // write out all the boundary types found
    ng.write_out_nodes(IBM_INNER, file_write);
    ng.write_out_nodes(IBM_OUTER, file_write);
    ng.write_out_nodes(NO_BOUNDARY, file_write);
    // write out the markers
    ng.markers->write_out_markers(file_write);
    std::cout << ng.markers->marker_points.size() << std::endl;
    // write out surface
    sg.write_out_surface();
    // end
    return 0;
}
