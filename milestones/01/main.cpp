#include <iostream>
#include "lattice_boltzmann.h"
#include <ctime>
#ifdef USE_MPI
#include <mpi.h>
#endif

int main(int argc, char *argv[]) {
    /// init + run functional test
    /// serach algorithm is hot garbage
    int size_x = 100;
    int size_y = 100;
    simulation sim;
    auto time = clock();
    sim.init(size_x, size_y);
    sim.run();
    time = clock() -time;
    std::cout << double(time)/CLOCKS_PER_SEC << std::endl;
}

