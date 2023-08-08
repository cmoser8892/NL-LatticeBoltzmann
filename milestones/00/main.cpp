
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
    long canvas_size = 100;
    kernelType_t kernel = KERNEL_C;
    double ibm_distance = kernel_id_to_lattice_search(kernel);
    // Load the image
    auto test_image = get_base_path();
    test_image.append("tests");
    test_image.append("test_images");
    test_image.append("test_small.png");
    // call the drawer
    surfaceDrawer s(test_image);
    // std::vector<int> sel = {0,4,8,12};
    std::vector<int> sel = {0};
    s.run_selective(sel);
    s.surface_storage.surface_mass_center();
    nodeGenerator ng(&s.surface_storage);
    ng.init_surface_return(canvas_size,ibm_distance);
    // write out all the boundary types found
    ng.write_out_nodes(IBM_INNER, true);
    ng.write_out_nodes(IBM_OUTER, true);
    ng.write_out_nodes(NO_BOUNDARY, true);
    // end
    return 0;
}
