#include <iostream>
#include "lattice_boltzmann.h"
#include <ctime>
#ifdef USE_MPI
#include <mpi.h>
#endif

int main(int argc, char *argv[]) {
    /// init + run functional test for timing and stuff
    int size_x = 3;
    int size_y = 3;
    simulation sim;
    auto time = clock();
    sim.init(size_x, size_y);

    std::cout << double(clock()-time)/CLOCKS_PER_SEC << std::endl;
    time = clock();
    for(int i = 0; i < 1000; ++i)
        sim.run(); // takes about a min
    std::cout << double(clock()-time)/CLOCKS_PER_SEC << std::endl;
    // prob takes way to long in the way that it is programmed right now
    // mainly memory bound
    sim.get_data();
}

