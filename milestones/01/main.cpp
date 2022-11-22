#include <Eigen/Dense>
#include <iostream>
#include "lattice_boltzmann.h"
#ifdef USE_MPI
#include <mpi.h>
#endif


int main(int argc, char *argv[]) {
    int size_x = 5;
    int size_y = 5;
    simulation sim;
    sim.init(size_x,size_y);
    for( int i = 0; i < 5; ++i) {
        sim.run();
    }
}
