#include "simulation.h"
#include <ctime>
#include <iostream>

int main(int argc, char *argv[]) {
    int size = 102;
    int steps = 100000;
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
    }
    sim.get_data(true);
    return 0;
}
