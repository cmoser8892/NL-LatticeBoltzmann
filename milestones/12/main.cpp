#include <iostream>
#include <chrono>
#include "lbm_simulation.h"
#include "drawn_image_surface.h"

/**
 * 12 main. We can change an image to be the basis of a simulation with ibm.
 * @attention Because ibm simulations appear to be unstable (could also be that i dont relax them long enough) only circular forces should be considered for now.
 * @note to change between used forces go the the ibmClass, you have to change it directly.
 * @param argc
 * @param argv
 * @return
 */
int main(int argc, char *argv[]) {
    int steps = 10000;
    long canvas_size = 200;
    double marker_distance = 0.5;
    kernelType_t kernel = KERNEL_C;
    // Load the image
    auto test_image = get_base_path();
    test_image.append("tests");
    test_image.append("test_images");
    test_image.append("donut.png");
    // image drawer stuff
    surfaceDrawer drawer(test_image);
    std::vector<int> sel = {0,4,8,12};
    drawer.run_selective(sel);
    drawer.surface_storage.surface_mass_center();
    nodeGenerator ng(&drawer.surface_storage);
    double ibm_distance = kernel_id_to_lattice_search(kernel);
    /// todo main diff between 11 and 11b kept in 12
    ng.init_surface_return(canvas_size,ibm_distance,marker_distance);
    ng.visualize_2D_nodes();
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
                    std::cerr << m->position.x() << ", " << m->position.y() << "//"
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
