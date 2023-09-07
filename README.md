# Code for the Masterthesis: Arbitary shapes in Lattice Boltzmann

The main goal is to create a datastructure that enables flexibilty with boundaries in a Lattice Boltzmann Simulation. 
Neighborhood lists with singular nodes are used instead of the usual grids to enable simulations with just the requiered amount of fluid nodes.
This work here explore this idea for structures with solid/liquid << 1. 
It applies opencv and other more manual methods to get from a simple drawing (can be any kind of image) to a simulation more or less seamlessly.
Supported are bounce back boundaries and IBM boundaries.

This was modified from the [HPC MD
with C++
project](https://imtek-simulation.github.io/MolecularDynamics/_project/general_remarks.html) made 
with 
 [CMake](https://cmake.org/).
 That github can be found here [HPC MD Github sekelton](https://github.com/IMTEK-Simulation/cmake-skeleton).
