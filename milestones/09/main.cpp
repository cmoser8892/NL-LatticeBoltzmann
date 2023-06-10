#include "simulation.h"
#include "image_converter.h"
#include <iostream>
#include <chrono>

int main(int argc, char *argv[]) {
    auto start = std::chrono::high_resolution_clock::now();
    // init variant
    int steps = 20000;
    // test image setupa
    auto bmp_24_test_image = get_base_path();
    bmp_24_test_image.append("tests");
    bmp_24_test_image.append("test_images");
    bmp_24_test_image.append("small_t.bmp");
    imageConverter ic(bmp_24_test_image);
    // run the functions
    ic.init();
    ic.run();
    // ic.boundaries->visualize_2D_boundary();
    // ic.boundaries->visualize_2D_boundary();
    nodeGenerator gen(ic.boundaries);
    // if the fused init runs this test is considered complete
    gen.init_fused(ic.return_basic_size());
    // init sim parameters
    simulation_parameters params;
    params.relaxation = 0.5;
    point_t dk = {0,0};
    goaForce rot(dk,ic.boundaries->size,6e-3);
    forcedSimulation sim(ic.boundaries,&gen, &rot);
    sim.set_simulation_parameters(params);
    sim.init();
    // run sim
    for(int i = 0; i < steps; ++i) {
        if(i % 1000 == 0) {
            std::cout << "Step: " << i << std::endl;
        }
        sim.run(i);
    }
    // write out the data
    sim.get_data(true);
    if(gen.straight_surfaces != nullptr)
        gen.straight_surfaces->write_out_surface();
    // time the whole thing
    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::seconds>(stop - start);
    std::cout << "Took " <<duration.count()<< "s" << std::endl;
    return 0;
}