
/*
 * vallgrind call
valgrind --tool=callgrind --dump-instr=yes (p)
 */

#include "helper_functions.h"
#include "lbm_simulation.h"
#include "drawn_image_surface.h"
#include "nodeGenerator.h"

/**
 * 00 main, testing stuff.
 * @note good rule of thumb is all 4 contours need to be selected
 * @return
 */
int main() {
    // node generator variables
    long canvas_size = 400;
    vector_t draw_size = {canvas_size-1,canvas_size-1};
    double marker_distance = 0.5;
    bool file_write = true;
    kernelType_t kernel = KERNEL_C;
    double ibm_distance = kernel_id_to_lattice_search(kernel);
    // Load the image
    auto test_image = get_base_path();
    test_image.append("tests");
    test_image.append("test_images");
    test_image.append("channel.png");
    // call the drawer
    surfaceDrawer s(test_image);
    std::vector<int> sel = {0};
    s.run_non_connecting(sel, true);
    // s.close_open_surface(draw_size);
    s.surface_storage.surface_mass_center();
    nodeGenerator ng(&s.surface_storage);
    ng.init_surface_return(canvas_size,ibm_distance,marker_distance);
    // new part
    // add additional surface
    if(0) {
        point_t p = {94,13};
        straight_t surface;
        surface.point = p;
        surface.direction = {0,1};
        surface.max_t = 32;
        surface.type = FAKE_FORCEING;
        // s.surface_storage.add_surface(surface);
        // detect and put markers on
        markerDistribution force;
        force.individual_distribute_markers(surface);
        // reflag nodes
        ng.reflag_force_nodes(&force,kernel_id_to_lattice_search(kernel));
    }
    // write out all the boundary types found
    ng.write_out_nodes(IBM_INNER, file_write);
    ng.write_out_nodes(IBM_OUTER, file_write);
    ng.write_out_nodes(NO_BOUNDARY, file_write);
    ng.write_out_nodes(FAKE_FORCEING,file_write);
    ng.write_out_nodes(FAKE_FORCEING_INNER,file_write);
    // write out the markers
    ng.markers->write_out_markers(file_write);
    s.surface_storage.write_out_surface();
    // end init basics now add to sim

    if(0) {
        // start sim
        simulation_parameters params;
        int steps = 30000;
        params.relaxation = 0.5;
        params.ibm_range = kernel_id_to_lattice_search(kernel);
        params.kernel_in_use = kernel;
        params.k = 1;
        point_t dk = {0,0};
        // max rotation is 7.5e-3
        vector_t sizes = {canvas_size,canvas_size};
        goaForce rot(dk,sizes,1e-3); // i dont really need the force
        ibmSimulation sim(&ng, &rot,ng.markers,sizes);
        // init
        sim.set_simulation_parameters(params);
        sim.init();
        // sim.add_force_markers(&force,0.001);
        for(int i = 0; i < steps; ++i) {
            if (i % 100 == 0) {
                std::cout << "Step: " << i << std::endl;
            }
            sim.run(i);
        }
        sim.get_data(true);
    }

    // end
    return 0;
}
