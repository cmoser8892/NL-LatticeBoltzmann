//
// Created by christoph on 24.01.23.
//

#include "neighborhood.h"
#include <iostream>

void orderingNodes::order(std::vector<nodePoint_t*> nodes) {
   // compute grid cell coordinates and store indices
   std::vector<coordinate_t> coordinates;
   std::vector<uint64_t> store_indices; // full 64 bits aka 8 bytes
   for(auto node  : nodes) {
       // we round down to make little cells of full integers
       coordinate_t coordinate;
       coordinate.x = std::floor(node->position.x());
       coordinate.y = std::floor(node->position.y());
       // std::cout << coordinate.x << " "<< coordinate.y << std::endl;
       coordinates.push_back(coordinate);
       // we will mask out some bits so better make sure that is not bigger than max
       assert(coordinate.x >= pow(2,22));
       assert(coordinate.y >= pow(2,22));
       // interleave the bits
   }
}