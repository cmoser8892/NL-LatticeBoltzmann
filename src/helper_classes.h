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

// watchdog for rho
class rhoWatchdog {
    // just prints out std error msgs
  private:
    flowfield_t rho;
    double sensitivity = 0.1; // displacement to the previous value
  public:
    rhoWatchdog(double s,point_t size);
    bool check(node* n,int step);
};

// works with bit interleaving floors points to a full point on the gird!
class pointKeyHash {
  private:
    std::unordered_multimap<handle_t,handle_t> keys;
  public:
    void fill_key(handle_t positions_handle, point_t pos);
    void clear();
    handle_t key_translation(point_t pos);
    handle_t key_translation(coordinate_t cord);
};

// basic stash that holds items for a number of iterations
// a window with a conveyor belt behind it, specifically for handles
class windowedHandles{
    // dont make this to big
  private:
    unsigned long target_size;
    std::list<handle_t> previous;
  public:
    // functions
    windowedHandles(unsigned long size);
    unsigned long size();
    void add(handle_t h);
    bool check(handle_t h);
};

#endif // NL_LATTICEBOLTZMANN_HELPER_CLASS_H
