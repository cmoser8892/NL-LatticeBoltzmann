//
// Created by christoph on 24.01.23.
//

#ifndef NL_LATTICEBOLTZMANN_NEIGHBORHOOD_H
#define NL_LATTICEBOLTZMANN_NEIGHBORHOOD_H
/*
now to actual neighbourhood stuff:
3 ideas:
- sort via cell indexes (like in the md class)
- z-shape ordering (used in particle simulations, and my favorite)
- hashed mapping
 --> hard question which algorithm has the best query time to the neighbors
 */

#include "types.h"
#include "node.h"
/**
 * first we put nodes into a cell
 * this may seem unnecessary for computation of an uniform grid
 * but may safe the day later if the node coordinates are not perfectly on int positions
 * for subdividing
 */
// order nodes into the cells
class neighbourhood {
  private:
    std::unordered_multimap<handle_t,handle_t> keys;
    void fill_keys(std::vector<nodePoint_t*> &nodes);
  public:
    void determine_neighbors(std::vector<nodePoint_t*> &nodes);
};


#endif // NL_LATTICEBOLTZMANN_NEIGHBORHOOD_H
