
/*
 * vallgrind call
valgrind --tool=callgrind --dump-instr=yes (p)
 */

#include "helper_functions.h"
#include "drawn_image_surface.h"
#include "nodeGenerator.h"

/*
 * todo implement the working plan:
 * Working plan here is to place a straight line into the surfaces
 * and cut every surface following it, they should be in the right order
 * till the next point is found
 * need to implement the rhoin+out here to but that should be a nonfactor
 */
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
    // Load the image
    auto test_image = get_base_path();
    test_image.append("tests");
    test_image.append("test_images");
    test_image.append("periodic_introduction.png");
    // call the drawer
    surfaceDrawer s(test_image);
    std::vector<int> sel = {0};
    s.run_non_connecting(sel, true);
    // s.close_open_surface(draw_size);
    // add additional inlet + outlet
    if(1) {
        straight_t inlet;
        inlet.type = PERIODIC;
        inlet.point = {165,160};
        inlet.direction = {0,1};
        inlet.max_t = 25;
        s.surface_storage.add_surface(inlet);
        straight_t outlet;
        outlet.type = PERIODIC;
        outlet.point = {93,115};
        outlet.direction = {1,0};
        outlet.max_t = 20;
        s.surface_storage.add_surface(outlet);
    }
    s.surface_storage.periodic_check_in();
    s.surface_storage.surface_mass_center();
    nodeGenerator ng(&s.surface_storage);
    ng.init_surface_return(canvas_size,kernel,marker_distance);
    // write out all the boundary types found
    ng.write_out_nodes(IBM_INNER, file_write);
    ng.write_out_nodes(IBM_OUTER, file_write);
    ng.write_out_nodes(NO_BOUNDARY, file_write);
    // write out the markers
    // ng.markers->write_out_markers(file_write);
    // std::cout << ng.markers->marker_points.size() << std::endl;
    s.surface_storage.write_out_surface();
    // end
    return 0;
}
