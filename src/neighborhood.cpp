#include "neighborhood.h"
#include "helper_functions.h"

/**
 * @fn neighbourhood::neighbourhood()
 * @brief constructor, sets up min and max coordinate
 */
neighbourhood::neighbourhood() {
    /// set the min to max and max to min
    min_coordinate.x = INT32_MAX;
    min_coordinate.y = INT32_MAX;
    max_coordinate.x = 0;
    max_coordinate.y = 0;
}

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
       // snoop the lower positions
       snoop_min_coordinate(coordinate);
       snoop_max_coordinate(coordinate);
       // fill the pkh
       pkh.fill_key(node->handle,node->position);
   }
}

/**
 * @fn void neighbourhood::snoop_min_coordinate(coordinate_t coordinate)
 * @brief snoops the min coorodinates and sets that accoridnglly
 * @param coordinate
 */
void neighbourhood::snoop_min_coordinate(coordinate_t coordinate) {
    // set the min coordinate to be reshuffled to
    if(coordinate.x < min_coordinate.x) {
        min_coordinate.x = coordinate.x;
    }
    if(coordinate.y < min_coordinate.y ) {
        min_coordinate.y = coordinate.y;
    }
}

/**
 * @fn void neighbourhood::snoop_max_coordinate(coordinate_t coordinate)
 * @brief snoops the max coordiante and sets that
 * @param coordinate
 */
void neighbourhood::snoop_max_coordinate(coordinate_t coordinate) {
    // set the max coordinate to be reshuffled to
    if(coordinate.x > max_coordinate.x) {
        max_coordinate.x = coordinate.x;
    }
    if(coordinate.y > max_coordinate.y ) {
        max_coordinate.y = coordinate.y;
    }
}

/**
 * @fn void neighbourhood::periodic_coordinate_reshuffle(coordinate_t* coordinate)
 * @brief used to setup links for periodic wet boundaries
 * @param coordinate
 */
void neighbourhood::periodic_coordinate_reshuffle(coordinate_t* coordinate) {
    // todo write a test to check the reshuffle
    // min set
    if(coordinate->x < min_coordinate.x) {
        // set to max
        coordinate->x = max_coordinate.x;
    }
    if(coordinate->y < min_coordinate.y) {
        // set to max
        coordinate->y = max_coordinate.y;
    }
    // max set
    if(coordinate->x > max_coordinate.x) {
        // set to min
        coordinate->x = min_coordinate.x;
    }
    if(coordinate->y > max_coordinate.y) {
        // set to min
        coordinate->y = min_coordinate.y;
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
           // extra coordinate handling
           if((node->boundary == PERIODIC) || (node->boundary == PRESSURE_PERIODIC)) {
                periodic_coordinate_reshuffle(&coordinate);
                current = {coordinate.x,coordinate.y};
           }
           handle_t found_handle = pkh.key_translation(coordinate);
           if(found_handle > 0) {
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
                       // add if a wet node
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

/**
 * @fn void neighbourhood::check_wet_nodes(std::vector<nodePoint_t *> &nodes)
 * @brief checks wet nodes channels, aka are all boundaries there
 * @param nodes
 */
void neighbourhood::check_wet_nodes(std::vector<nodePoint_t *> &nodes) {
    // checks wet link size should be Channels-1 so 8 for D2Q9
    for(auto node: nodes) {
        if(node->type == WET) {
            if(node->links.size() != CHANNELS-1) {
                throw std::invalid_argument("Invalid Channels");
            }
        }
    }
}