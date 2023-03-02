#ifndef NL_LATTICEBOLTZMANN_NEIGHBORHOOD_H
#define NL_LATTICEBOLTZMANN_NEIGHBORHOOD_H

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
