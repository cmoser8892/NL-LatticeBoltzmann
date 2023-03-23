#ifndef NL_LATTICEBOLTZMANN_NEIGHBORHOOD_H
#define NL_LATTICEBOLTZMANN_NEIGHBORHOOD_H

#include "types.h"
#include "node.h"
#include "helper_functions.h"
/**
 * first we put nodes into a cell
 * this may seem unnecessary for computation of an uniform grid
 * but may safe the day later if the node coordinates are not perfectly on int positions
 * for subdividing
 */
// order nodes into the cells
class neighbourhood {
  private:
    pointKeyHash pkh;
    coordinate_t min_coordinate;
    coordinate_t max_coordinate;
    void fill_keys(std::vector<nodePoint_t*> &nodes);
    void snoop_min_coordinate(coordinate_t coordinate);
    void snoop_max_coordinate(coordinate_t coordinate);
    void periodic_coordinate_reshuffle(coordinate_t* coordinate);
  public:
    neighbourhood();
    void determine_neighbors(std::vector<nodePoint_t*> &nodes);
    void check_wet_nodes(std::vector<nodePoint_t*> &nodes);
};


#endif // NL_LATTICEBOLTZMANN_NEIGHBORHOOD_H
