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
typedef struct orderedCell {
    unsigned x;
    unsigned y;
    std::vector<handle_t> inhabitants;
}orderedCell_t;

// order nodes into the cells
class orderingNodes {
  private:
    std::vector<orderedCell_t*> cells;
  public:
    void order(std::vector<nodePoint_t*> nodes);
};


#endif // NL_LATTICEBOLTZMANN_NEIGHBORHOOD_H
