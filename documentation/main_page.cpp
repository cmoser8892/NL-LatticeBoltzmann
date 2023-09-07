// note this is just cpp cause it makes clion help me write, not part of the code!!

/**
 * \mainpage Lattice Boltzmann with Neighborhood Lists
 * \subsection General
 * The project is written in C++20. The main goal is to create a datastructures
 * that enables us to be quite flexible with boundaries in a Lattice Boltzmann Simulation.
 * For this we use neighborhood lists to connect a bunch of loose populations together.
 * This work explores that idea.
 * \subsection Dependencies
 * Eigen 3.10. \n
 * C++ Standard Libraries \n
 * Doxygen version 1.8.17 (It looks different if you use a different version duh)
 * \subsection Limitations
 * The Neighborhood list sorts every node into its own cell,
 * great if you want a simple algorithm but bad for high node counts. \n
 * Only positive coordinates will work for any kind of structure.
 * (The hashes in the nl must not be negative (64 bit hashes, constructed from the interleaved floored position, negative numbers are unaccounted for)).\n
 * (This results in a maximum canvas draw size of 0 to 2e32(4gb) in 2d and 0 to 2e21(2mb) in 3d)
 * Handles to the nodes start at 1.
 * \subsection Ideas
 * Lattice Boltzmann with Neighborhood lists seems cool for strongly moving boundaries.
 *  It would be a very easy to implement a grassroot way from the nodes up instead of the hodgepotch i have now (first try is messy).
 *  (Mix between nodes and markers). We can also to the inefficient way of going from the nodes to the markers
 *  and not the other way around, making everything fully node controlled and based on the ns algorithm.
 *  Using a pre ibm node as a marker and checking weather or not one of the pre markers is actually a valid population that has to created. \n
 *  Would be also cool to try out an extreme ns based approach, right now i build the list and use a static link list in the nodes
 *  \subsection Particularities
 *  Tests are messy and not really organized in files but in the Test suits (that is why there is a number of misc-tests, I only partition them when the file gets too big).
 *  Also i would describe them as a bit of a italian mess, there is a lot of copy paste, but in my opinion a test should always be self complete.
 *  Everything undocumented in doxygen should be regarded with suspicion. In my normal workflow I document a method after testing.
 *  Anything undocumented can therefore be regarded as at least untested, maybe even non-functional.
 *  Critique on programming is very welcome :)
*/

