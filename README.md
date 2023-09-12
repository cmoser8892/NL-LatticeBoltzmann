
# Code to the Masterthesis - Lattice Boltzmann Method: Arbitary shapes

The main goal is to create a datastructure that enables flexibilty with boundaries in a Lattice Boltzmann Simulation. 
Neighborhood lists with singular nodes are used instead of the usual grids to enable simulations with just the requiered amount of fluid nodes.
This work here explore this idea for structures with solid/liquid << 1. 
It applies OpenCV and other more manual methods to get from a simple drawing (can be any kind of image) to a simulation more or less seamlessly.
Supported are bounce back boundaries and IBM boundaries.

The experineced reader my dircly skip into the documentation folder of the project and generate the doxygen to the project. Here a brief overview over each folder is given to orientate anyone.

| Folder: | Description: |
|----------|----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| documentation |  The generation of the Doxygen documentation is described here. After being generated the files for the documentation are also located in this folder. |
| milestones | The individual simulations shown in chapter \ref{ch:exp}, going from sliding lids over to variable structures simulated are located.  Sub variants are denoted with letters, for example 11a, which is a variant of milestone 11. For individual problems tackled it is best to check the Doxygen documentation. |
| simulation |  This folder contains the continuous development of the simulation class. To ensure backwards compatibility once a specific simulation would have been to dissimilar to previous work in the milestones a new class was created. This was also done to reduce the bloat in this class as the in initial idea rarely was the best and was incremented upon. They represent the development form a more basic simulation with a two-step algorithm to the optimized one-step algorithm. Next forces were introduced and a new class was created. To limit the complexity and just use the best parts, it was once again rewritten into lbm-simulations. |
| src | In this folder all the other source code is located. It includes all the classes necessary for the initialization with lots of helper function and classes used to initialize the relevant nodes. For a reference and explanations it is best to look into the code directly. |
| tests | All the tests are ordered within test suites and no one file includes all the tests of each test suite. As a reason the author refers to the fact that any file with more than 1000 lines tended to be slow to load. In general each test should be self contained with no additional functions needed but the ones for setup and testing. The individual test suits are: FunctionalTest where the simulation class, helper functions and classes are tested. All the force interactions, initialization in general, streaming, tests related to bounce back and the testing of IBM boundaries have their own suit too. |


The seklton used was modified from the [HPC MD
with C++
project](https://imtek-simulation.github.io/MolecularDynamics/_project/general_remarks.html) made 
with 
 [CMake](https://cmake.org/).
 That github can be found here [HPC MD Github sekelton](https://github.com/IMTEK-Simulation/cmake-skeleton).
