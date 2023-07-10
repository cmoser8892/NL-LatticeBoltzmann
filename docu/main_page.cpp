// note this is just cpp cause it makes clion help me write, not part of the code!!

/** \mainpage Lattice Boltzmann with Neighborhood Lists
 * \subsection General
 * The project is written in C++20. The main goal is to create a datastructures
 * that enables us to be quite flexible with boundaries in a Lattice Boltzmann Simulation.
 * For this we use neighborhood lists to connect a bunch of loose populations together.
 * \subsection Dependencies
 * Eigen 3.10. \n
 * C++ Standard Libraries \n
 * Doxygen version 1.8.17
 * \subsection Limitations
 * The Neighborhood list sorts every node into its own cell,
 * great if you want a simple algorithm but bad for high node counts. \n
 * Only positive coordinates will work for any kind of structure.
 * (The hashes in the nl must not be negative (64 bit hashes, constructed from the interleaved floored position, negative numbers are unaccounted for)).\n
 * Handles are start at 1.
 *
*/

