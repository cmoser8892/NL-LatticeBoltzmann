#ifndef NL_LATTICEBOLTZMANN_HELPER_CLASS_H
#define NL_LATTICEBOLTZMANN_HELPER_CLASS_H

#include "types.h"
#include "node.h"
#include "helper_functions.h"
#include "functions.h"
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
    bool check_force(fNode* n, int step);
};

/**
 * Watchdog to make sure our markers do not move too much, will yell when they do.
 */
class markerWatchdog {
  private :
    std::vector<point_t> previous; /**< container of all the previous positions */
    std::vector<point_t> original; /**< container of all the original positions */
    double sensitivity = 0.1; /**< sensitivity 10% */
  public:
    explicit markerWatchdog(double s);
    void init_marker(point_t p);
    bool check(point_t p, handle_t pos);
};

/**
 * Hashes points to be found in a neighborhood list.
 * @note works with bit interleaving floors points to a full point on the gird!
 * @attention For bigger problems this is the hard bottleneck in the init, cells of size 1 are not ideal here.
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
 * Boilerplate toe hide some of the ranging functionality of ranging points.
 * Not like the pgk this class only returns the valid point handles that are in range, that we search in.
 * @note not part of the pointKeyHash by choice, we do not want that in there
 * @todo gradual support for all the different type of nodes/points for false save mode
 * @attention main difference to pkh is that there is a guaranty that only handles in the range get returned
 */
class rangingPointKeyHash {
  private:
    pointKeyHash pkh; /**< Point key hash stash */
    std::vector<point_t> points;  /**< Points in the stash*/
    bool safe_position_yes_no = true; /**< should be set to true if we want to safe the points */
  public:
    void set_position_save(bool set);
    void fill_key(handle_t position_handle, point_t position);
    void clear();
    long size();
    std::vector<handle_t> ranging_key_translation(point_t pos, double range);
    bool ranging_key_look_for_specific(point_t pos, double range, handle_t h);
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
