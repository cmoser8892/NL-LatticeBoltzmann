#include "neighborhood.h"
#include "helper_functions.h"

/**
 * @fn void neighbourhood::fill_keys(std::vector<nodePoint_t*> &nodes)
 * @brief fills the keys map based on the z-order of the position making a hash map
 * @param nodes
 */
void neighbourhood::fill_keys(std::vector<nodePoint_t*> &nodes) {
   // compute grid cell coordinates and store keys
   for(auto node  : nodes) {
       // we round down to make little cells of full integers
       coordinate_t coordinate;
       coordinate.x = std::floor(node->position.x());
       coordinate.y = std::floor(node->position.y());
       // we will mask out some bits so better make sure that is not bigger than max
       assert(coordinate.x <= pow(2,22));
       assert(coordinate.y <= pow(2,22));
       // interleave the bits and write it into the table
       handle_t key = bit_interleaving_2d(coordinate.x,coordinate.y);
       keys.emplace(key,node->handle);
   }
}

/**
 * @fn void neighbourhood::determine_neighbors(std::vector<nodePoint_t *> &nodes)
 * @brief calculates the z-order of the position of the neighbours and looks for them in the hash table
 * @param nodes
 */
void neighbourhood::determine_neighbors(std::vector<nodePoint_t *> &nodes) {
   // fill the keys
   fill_keys(nodes);
   for(auto node: nodes) {
       // loop through the channels
       for(int i = 1; i < CHANNELS; ++i) {
           // set the position the the corressponding coordinate
           point_t current = node->position + velocity_set.col(i);
           coordinate_t coordinate;
           coordinate.x = std::floor(current.x());
           coordinate.y = std::floor(current.y());
           // interleave bits to get the correct key for that nodes handle
           handle_t search_key = bit_interleaving_2d(coordinate.x,coordinate.y);
           // try to find the key
           if(auto found_iter = keys.find(search_key); found_iter != keys.end()) {
               // set handle and array position remember handles start at 1!
               handle_t found_handle = found_iter->second;
               handle_t array_position = found_handle -1;
               // setup temp test point
               auto temp = point_t(nodes.at(array_position)->position);
               // check if the found nodes key variables match
               if (compare_two_points(&current, &temp)) {
                   // always include all neighbours if we are a wet node
                   bool add_me = false;
                   if (node->type == WET) {
                       add_me = true;
                   }
                   // if we are a dry node only include nodes that are wet
                   else if (node->type == DRY) {
                       if (nodes.at(array_position)->type == WET) {
                           add_me = true;
                       }
                   }
                   // if one of the above conditions holds add
                   if (add_me) {
                       toLinks_t link;
                       link.handle = found_handle;
                       link.channel = i;
                       node->links.push_back(link);
                   }
               }
           }
       }
   }
}

void neighbourhood::check_wet_nodes(std::vector<nodePoint_t *> &nodes) {
    // checks wet link size should be Channels-1 so 8 for D2Q9
    // todo no everyone has 8 channels special cases in the tests
    for(auto node: nodes) {
        if(node->type == WET) {
            if(node->links.size() != CHANNELS-1) {
                throw std::invalid_argument("Invalid Channels");
            }
        }
    }
}