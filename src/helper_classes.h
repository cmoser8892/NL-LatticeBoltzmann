#ifndef NL_LATTICEBOLTZMANN_HELPER_CLASS_H
#define NL_LATTICEBOLTZMANN_HELPER_CLASS_H

#include "types.h"
#include "node.h"
#include "helper_functions.h"
#include <iostream>
#include <sstream>
#include <vector>
#include <list>
#include <filesystem>

/**
 * Rho watchdog that makes sure rho values do not run away.
 */
class rhoWatchdog {
    // just prints out std error msgs
  private:
    flowfield_t rho; /**< Previous values */
    double sensitivity = 0.1; /**<  maximum displacement to the previous value */
  public:
    rhoWatchdog(double s,point_t size);
    bool check(node* n,int step);
};

/**
 * Hashes points to be found in a neighborhood list.
 * @note works with bit interleaving floors points to a full point on the gird!
 * @details works by reducing the position to a normal number and interleaving those number to create a classifier for that handle
 */
class pointKeyHash {
  private:
    std::unordered_multimap<handle_t,handle_t> keys; /**< handle handle parings */
   public:
    void fill_key(handle_t positions_handle, point_t pos);
    void clear();
    long map_size();
    handle_t key_translation(point_t pos);
    handle_t key_translation(coordinate_t cord);
    std::vector<handle_t> multi_key_translation(point_t pos);
    std::vector<handle_t> multi_key_translation(coordinate_t cord);
    std::vector<handle_t> ranging_key_translation(point_t pos, double range);
    std::vector<handle_t> ranging_key_translation(coordinate_t coord, double range);
};

/**
 * Basic stash that holds items for a number of iterations.
 * A window with a conveyor belt behind it, specifically for handles.
*/
class windowedHandles{
    // dont make this to big
  private:
    unsigned long target_size; /**<  size of the window that holds the handles */
    std::list<handle_t> previous; /**< vector that holds that past handles */
  public:
    // functions
    windowedHandles(unsigned long size);
    unsigned long size();
    void add(handle_t h);
    bool check(handle_t h);
};

#endif // NL_LATTICEBOLTZMANN_HELPER_CLASS_H
