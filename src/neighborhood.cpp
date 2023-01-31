//
// Created by christoph on 24.01.23.
//

#include "neighborhood.h"
#include "helper_functions.h"
#include <iostream>
#include <algorithm>

/// todo curb down
// order nodes base on the z spacial curve
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
       store_indices.push_back(bit_interleaving_2d(coordinate.x,coordinate.y));
   }
   assert(coordinates.size()==store_indices.size());
   assert(coordinates.size() == nodes.size());
   // rewrite the handles
   for(int i = 0; i < nodes.size();++i) {
       // we write z order + 1 as a new handle
       nodes.at(i)->handle = store_indices.at(i) + 1;
   }
   // sort
   std::sort(nodes.begin(),nodes.end(), compare_handles);
}