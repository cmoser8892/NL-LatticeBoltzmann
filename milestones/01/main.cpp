#include <iostream>
#include "lattice_boltzmann.h"
#include <ctime>
#ifdef USE_MPI
#include <mpi.h>
#endif

int main(int argc, char *argv[]) {
    /// init + run functional test
    /// serach algorithm is hot garbage
    int size_x = 300;
    int size_y = 300;
    simulation sim;
    auto time = clock();
    sim.init(size_x, size_y);
    std::cout << double(clock()-time)/CLOCKS_PER_SEC << std::endl;
    sim.run();
    std::cout << double(clock()-time)/CLOCKS_PER_SEC << std::endl;
}

