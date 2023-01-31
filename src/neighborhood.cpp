//
// Created by christoph on 24.01.23.
//

#include "neighborhood.h"
#include "helper_functions.h"
#include <iostream>
#include <algorithm>

/// todo curb down
// order nodes base on the z spacial curve
void orderingNodes::order(std::vector<nodePoint_t*> &nodes) {
   // compute grid cell coordinates and store indices
   for(auto node  : nodes) {
       // we round down to make little cells of full integers
       coordinate_t coordinate;
       coordinate.x = std::floor(node->position.x());
       coordinate.y = std::floor(node->position.y());
       // std::cout << coordinate.x << " "<< coordinate.y << std::endl;
       // we will mask out some bits so better make sure that is not bigger than max
       assert(coordinate.x <= pow(2,22));
       assert(coordinate.y <= pow(2,22));
       // interleave the bits and write it into the handle
       node->handle = bit_interleaving_2d(coordinate.x,coordinate.y);
   }
   // sort
   std::sort(nodes.begin(),nodes.end(), compare_handles);
   // rewrite into valid handles
   handle_t valid_handle = 1;
   for(auto node : nodes) {
       node->handle = valid_handle++;
   }
}