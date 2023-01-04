#include "simulation.h"
#include <ctime>
#include <iostream>

int main(int argc, char *argv[]) {
    int size = 52;
    int steps = 10;
    point_t p = {size,size};
    boundaryPointConstructor boundaries(p);
    boundaries.init_sliding_lid();
    //
    simulation sim(&boundaries);
    sim.init();
    for(int i = 0; i < steps; ++i) {
        sim.run();
        if(i % 1000 == 0) {
            std::cout << "Step: " << i << std::endl;
        }
        sim.get_data(false);
    }
    sim.get_data(true);
    return 0;
}
