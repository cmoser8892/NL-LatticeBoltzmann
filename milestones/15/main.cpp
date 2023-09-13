
/*
 * vallgrind call
valgrind --tool=callgrind --dump-instr=yes (p)
 */

#include "helper_functions.h"
#include "drawn_image_surface.h"
#include "nodeGenerator.h"
#include "lbm_simulation.h"

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
    long canvas_size = 800;
    vector_t draw_size = {canvas_size-1,canvas_size-1};
    double marker_distance = 0.5;
    bool file_write = true;
    kernelType_t kernel = KERNEL_C;
    // Load the image
    auto test_image = get_base_path();
    test_image.append("tests");
    test_image.append("test_images");
    test_image.append("mod_snake_structure.png");
    // call the drawer
    surfaceDrawer s(test_image);
    std::vector<int> sel = {0};
    s.run_non_connecting(sel, true);
    // s.close_open_surface(draw_size);
    // add additional inlet + outlet
    if(1) {
        straight_t inlet;
        inlet.type = PERIODIC;
        inlet.point = {400,465};
        inlet.direction = {1,0};
        inlet.max_t = 100;
        s.surface_storage.add_surface(inlet);
        straight_t outlet;
        outlet.type = PERIODIC;
        outlet.point = {615,650};
        outlet.direction = {0,1};
        outlet.max_t = 100;
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
    std::cout << ng.markers->marker_points.size() << std::endl;
    s.surface_storage.write_out_surface();
    // setup params for sim
    int steps = 20000;
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
        sim.run(i);
    }
    //
    sim.get_data(true);
    if(ng.straight_surfaces != nullptr)
        ng.straight_surfaces->write_out_surface();
    // end
    return 0;
}
