
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
    long canvas_size = 200;
    vector_t draw_size = {canvas_size-1,canvas_size-1};
    double marker_distance = 0.5;
    bool file_write = true;
    kernelType_t kernel = KERNEL_C;
    double ibm_distance = kernel_id_to_lattice_search(kernel);
    // Load the image
    auto test_image = get_base_path();
    test_image.append("tests");
    test_image.append("test_images");
    test_image.append("donut.png");
    // call the drawer
    surfaceDrawer s(test_image);
    std::vector<int> sel = {0,4};
    s.run_non_connecting(sel, true);
    // s.close_open_surface(draw_size);
    s.surface_storage.surface_mass_center();
    nodeGenerator ng(&s.surface_storage);
    ng.init_surface_return(canvas_size,ibm_distance,marker_distance);
    // new part
    // add additional surface
    point_t p = {94,13};
    straight_t surface;
    surface.point = p;
    surface.direction = {0,1};
    surface.max_t = 32;
    surface.type = FAKE_FORCEING;
    s.surface_storage.add_surface(surface);
    // detect and put markers on
    markerIBM force;
    force.individual_distribute_markers(surface);
    // reflag nodes
    ng.reflag_force_nodes(&force,2);
    //
    // write out all the boundary types found
    ng.write_out_nodes(IBM_INNER, file_write);
    ng.write_out_nodes(IBM_OUTER, file_write);
    ng.write_out_nodes(NO_BOUNDARY, file_write);
    ng.write_out_nodes(FAKE_FORCEING,file_write);
    ng.write_out_nodes(FAKE_FORCEING_INNER,file_write);
    // write out the markers
    ng.markers->write_out_markers(file_write);
    s.surface_storage.write_out_surface();
    // end
    return 0;
}
