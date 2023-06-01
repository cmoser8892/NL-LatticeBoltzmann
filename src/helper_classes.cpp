#include "helper_classes.h"

/**
 * @fn rhoWatchdog::rhoWatchdog(double s,point_t size
 * @brief constructor for the rho_watchdog
 * @param s
 * @param size
 */
rhoWatchdog::rhoWatchdog(double s,point_t size) :sensitivity(s) {
    rho.setOnes(long(size.x()),long(size.y()));
}

/**
 * @fn bool rhoWatchdog::check(node *n,int step)
 * @brief performs a watchdog check of the history of the rho value, aka compares it to the previous one
 * @param n
 * @param step
 * @return
 */
bool rhoWatchdog::check(node *n,int step) {
    double rho_old = rho(int(n->position(0)),int(n->position(1)));
    bool return_value = false;
    if((abs(rho_old-n->rho)) >= (abs(rho_old*sensitivity))) {
        std::cerr << "Rho-diviation at " << step << std::endl;
        std::cerr << "Position: " << n->position.x()
                  << " ," << n->position.y()
                  << std::endl;
        std::cerr << "Rho previous: " << rho_old << std::endl;
        std::cerr << "Rho now: " << n->rho << std::endl;
        std::cerr << std::endl;
        return_value = true;
    }
    rho(int(n->position(0)),int(n->position(1))) = n->rho;
    return return_value;
}

/**
 * @fn void pointKeyHash::fill_key(handle_t positions_handle, point_t pos)
 * @brief put a position in the key table, floors the double value
 * @param positions_handle
 * @param pos
 */
void pointKeyHash::fill_key(handle_t positions_handle, point_t pos) {
    coordinate_t coordinate;
    coordinate.x = std::floor(pos.x());
    coordinate.y = std::floor(pos.y());
    handle_t key = bit_interleaving_2d(coordinate.x,coordinate.y);
    if(coordinate.x >= pow(2,22)) {
        throw std::invalid_argument("Overflow while bit-interleaving");
    }
    if(coordinate.y >= pow(2,22)) {
        throw std::invalid_argument("Overflow while bit-interleaving");
    }
    if(coordinate.x < 0 || coordinate.y < 0) {
        // if there is negative number in the hash table we get two entries which we can not handle
        throw std::invalid_argument("PKH has to be modified to handle negative numbers");
    }
    keys.emplace(key,positions_handle);
}

/**
 * @fn handle_t pointKeyHash::key_translation(point_t pos)
 * @brief finds the key or does not in this case returns 0
 * @param pos a point
 * @return
 */
handle_t pointKeyHash::key_translation(point_t pos) {
    handle_t return_key = 0;
    // generate the search key
    coordinate_t coordinate;
    coordinate.x = std::floor(pos.x());
    coordinate.y = std::floor(pos.y());
    return key_translation(coordinate);
}

/**
 * @fn handle_t pointKeyHash::key_translation(coordinate_t coord)
 * @brief coordinate based key translation
 * @param coord
 * @return
 */
handle_t pointKeyHash::key_translation(coordinate_t coord) {
    handle_t return_key = 0;
    // todo no handling if there is more than one result found
    handle_t search_key = bit_interleaving_2d(coord.x,coord.y);
    if(auto found_iter = keys.find(search_key); found_iter != keys.end()) {
        // does not do the translation into an array position
        return_key = found_iter->second;
    }
    return return_key;
}

/**
 * @fn void pointKeyHash::clear()
 * @brief clears the pkh key table
 */
void pointKeyHash::clear() {
    keys.clear();
}

/**
 * @fn windowedHandles::windowedHandles(unsigned long s)
 * @brief constructor set up the list with the wished for size
 * @param s
 */
windowedHandles::windowedHandles(unsigned long s) {
    // set the size for the vector container
    target_size = s;
    // init the list to the target size
    for(int i = 0; i < target_size; ++i) {
        previous.push_back(0);
    }
}

/**
 * @fn unsigned long windowedHandles::size()
 * @brief returns the size of the list data-structure
 * @return size of the list data-structure
 */
unsigned long windowedHandles::size() {
    return previous.size();
}

/**
 * @fn void windowedHandles::add(handle_t h)
 * @brief adds a new element into the list pops the oldest one
 * @param h
 */
void windowedHandles::add(handle_t h) {
    // pop the oldest element and push in a new one
    previous.push_front(h);
    previous.pop_back();
}

/**
 * @fn bool windowedHandles::check(handle_t h)
 * @brief checks if the given handle h is in the buffer
 * @param h
 * @return returns true if the handle h is in the buffer
 */
bool windowedHandles::check(handle_t h) {
    // answers the yes no question if previously seen
    bool return_value = false;
    for (auto element : previous) {
        if (element == h) {
            return_value = true;
            break;
        }
    }
    // return me
    return return_value;
}