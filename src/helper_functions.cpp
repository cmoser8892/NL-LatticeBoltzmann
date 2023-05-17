#include <fstream>
#include <iostream>
#include <x86intrin.h> // x86 instructions
#include "helper_functions.h"

/**
 * @fn bool node_position_comparison(node* n, array_t* position)
 * @brief compares the nodes postion with an postion in an array
 * @param n
 * @param position
 * @return
 */
bool node_position_comparison(node* n, array_t* position) {
    array_t* node_position = &n->position;
    // check same size
    if(node_position->size() != position->size()) {
        return false;
    }
    // check same values
    for(int i = 0; i < node_position->size(); ++i) {
        if(node_position->operator()(i) == position->operator()(i)) {}
        else {
            return false;
        }
    }
    return true;
}

/**
 * @fn void write_flowfield_data(flowfield_t * field, std::string filename, bool write_to_file)
 * @brief writes the data form a flow-field into a file
 * @param field
 * @param filename
 * @param write_to_file
 */
void write_flowfield_data(flowfield_t * field, std::string filename, bool write_to_file) {
    std::ofstream out;
    out.open(filename);
    // preamble size infos
    // print out the rows first then to file
    for(int i = 0; i < field->rows(); ++i) {
        if(write_to_file) {
            out << field->row(i) << std::endl;
        }
        else {
            std::cout << field->row(i) << std::endl;
        }
    }
    out.close();
}

/**
 * @fn bool check_inside_limits_upper_lower(point_t* p, point_t* limit_lower, point_t* limit_upper)
 * @brief checks weather or not a point is inside the limit defined py two points
 * @param p
 * @param limit_lower
 * @param limit_upper
 * @return
 */
bool check_inside_limits_upper_lower(point_t* p, point_t* limit_lower, point_t* limit_upper) {
    bool return_value = true;
    if(p->x() < limit_lower->x() || p->x() > limit_upper->x()) {
        return_value = false;
    }
    if(p->y() < limit_lower->y() || p->y() > limit_upper->y()) {
        return_value = false;
    }
    return return_value;
}
 /**
  * @fn bool compare_two_points(point_t* p1, point_t* p2 )
  * @brief compares tow points
  * @param p1
  * @param p2
  * @return
  */
bool compare_two_points(point_t* p1, point_t* p2 ) {
    // double comparison is actually == standard
    bool return_value = true;
    if( p1->x() != p2->x()) {
        return_value = false;
    }
    if( p1->y() != p2->y()) {
        return_value = false;
    }
    return return_value;
}
/**
 * @fn std::vector<std::string> split_string (std::string s, std::string delimiter)
 * @brief splits a string
 * @param s
 * @param delimiter
 * @return string parts in a vector
 */
std::vector<std::string> split_string (std::string s, std::string delimiter) {
    size_t pos_start = 0, pos_end, delim_len = delimiter.length();
    std::string token;
    std::vector<std::string> res;

    while ((pos_end = s.find (delimiter, pos_start)) != std::string::npos) {
        token = s.substr (pos_start, pos_end - pos_start);
        pos_start = pos_end + delim_len;
        res.push_back (token);
    }

    res.push_back (s.substr (pos_start));
    return res;
}

/**
 * @fn uint64_t bit_interleaving(uint32_t x, uint32_t y)
 * @brief classical bit interleaving without x86
 * @param x
 * @param y
 * @return interleaved
 */
uint64_t bit_interleaving(uint32_t x, uint32_t y) {
    uint64_t z = 0;
    for(int i = 0; i < CHAR_BIT*sizeof(x); ++i)
        z|=((x&(1<<i))<<i)|((y&(1<<i))<<(i+1));
    return z;
}

/**
 * @fn uint64_t bit_interleaving_2d(uint32_t x, uint32_t y)
 * @brief _pdep bit interleaving, specific to x86 thou
 * @param x
 * @param y
 * @return interleaved bits
 */
uint64_t bit_interleaving_2d(uint32_t x, uint32_t y) {
    return _pdep_u64(x,0x5555555555555555) | _pdep_u64(y,0xaaaaaaaaaaaaaaaa);
}

/**
 * @fn uint64_t bit_interleaving_3d(uint32_t x, uint32_t y, uint32_t z)
 * @brief interleaves 32 bits using the x86 instruction directly
 * @param x only least significant 21 bits are interleaved
 * @param y only least significant 21 bits are interleaved
 * @param z only least significant 21 bits are interleaved
 * @return interleaved result
 */
uint64_t bit_interleaving_3d(uint32_t x, uint32_t y, uint32_t z) {
    // mask out the first part
    uint32_t mask = 0x0001FFFFF;
    x &= mask;
    y &= mask;
    z &= mask;
    uint64_t c =   _pdep_u64(x,0x9249249249249249)
                 | _pdep_u64(y,0x2492492492492492)
                 | _pdep_u64(z,0x4924924924924924);
    return c;
}

/**
 * @fn uint32_t bit_extraleaving_2d_x(uint64_t v)
 * @brief extraleaving function with a specific mask
 * @param v
 * @return
 */
uint32_t bit_extraleaving_2d_x(uint64_t v) {
    return _pext_u64(v,0x5555555555555555);
}

/**
 * @fn uint32_t bit_extraleaving_2d_y(uint64_t v)
 * @brief extraleaving function with a specific mask
 * @param v
 * @return
 */
uint32_t bit_extraleaving_2d_y(uint64_t v) {
    return _pext_u64(v,0xaaaaaaaaaaaaaaaa);
}

/**
 * @fn uint32_t bit_extraleaving_3d_x(uint64_t v)
 * @brief extraleaving function with a specific mask
 * @param v
 * @return
 */
uint32_t bit_extraleaving_3d_x(uint64_t v) {
    return _pext_u64(v,0x9249249249249249);
}

/**
 * @fn uint32_t bit_extraleaving_3d_y(uint64_t v)
 * @brief extraleaving function with a specific mask
 * @param v
 * @return
 */
uint32_t bit_extraleaving_3d_y(uint64_t v) {
    return _pext_u64(v,0x2492492492492492);
}

/**
 * @fn uint32_t bit_extraleaving_3d_z(uint64_t v)
 * @brief extraleaving function with a specific mask
 * @param v
 * @return
 */
uint32_t bit_extraleaving_3d_z(uint64_t v) {
    return _pext_u64(v,0x4924924924924924);
}

/**
 * @fn uint32_t reduce_32_2(uint32_t b)
 * @brief takes the most significant and least significant bits and reduces them  to 2 bits
 * @param b
 * @return 2 bit number
 */
uint32_t reduce_32_2(uint32_t b) {
    return _pext_u32(b,0x80000001);
}

/**
 * @fn std::filesystem::path get_base_path()
 * @brief gets the path to the NL-directory
 * @return the path to the NL-directory
 */
std::filesystem::path get_base_path() {
    auto executable_path = std::filesystem::current_path();
    std::filesystem::path return_path;
    auto search_string = "NL-LatticeBoltzmann";
    for(auto s : executable_path) {
        return_path.append(s.string());
        if(s == search_string) {
            break;
        }
    }
    return return_path;
}

/**
 * @fn std::vector<vector_t> circular_force_generation(int total_steps, int switch_time, double magnitude)
 * @brief
 * @param total_steps
 * @param switch_time
 * @param magnitude
 * @return
 */
std::vector<vector_t> circular_force_generation(int total_steps, int switch_time, double magnitude) {
    //
    std::vector<vector_t> return_vector;
    matrix_t helper = {{1,0,-1,0},
                       {0,1,0,-1}};
    int selector = 0;
    for(int i = 0; i < total_steps; ++i) {
        // switch selector up and down
        if((total_steps % switch_time) == 0) {
            selector = (++selector) % 3;
        }
        // we select on of the vectors and multiply with the magnitude
        vector_t current = helper.col(selector) * magnitude;
        // push back
        return_vector.push_back(current);
    }
    return return_vector;
}
/// helper and sub-classes
/**
 * @fn circularForce::circularForce(long switchtime, double mag)
 * @brief constructor of the circular force calculation
 * @param switchtime
 * @param mag
 */
circularForce::circularForce(long switchtime, double mag) {
    switch_time = switchtime;
    magnitude = mag;
    counter = 1;
    // setup x and y force
    x_force = {{1,0,-1,0}};
    y_force = {{0,1,0,-1}};
    x_force *= magnitude;
    y_force *= magnitude;
    // setup selectors
    selector_switchero();
}

/**
 * @fn
 * @brief
 */
void circularForce::selector_switchero() {
    if((counter % switch_time) == 0) {
        current_selector =  (++current_selector)%3;
    }
    if(((counter +1) % switch_time) == 0) {
        next_selector =  (++next_selector)%3;
    }
}

/**
 * @fn double circularForce::return_current_x()
 * @brief returns th current x force component
 * @return
 */
double circularForce::return_current_x() {
    return x_force(current_selector);
}

/**
 * @fn double circularForce::return_current_y()
 * @brief returns the current y force component
 * @return
 */
double circularForce::return_current_y() {
    return y_force(current_selector);
}

/**
 * @fn double circularForce::return_next_x()
 * @brief returns the next x force component
 * @return
 */
double circularForce::return_next_x() {
    return x_force(next_selector);
}

/**
 * @fn double circularForce::return_next_y()
 * @brief returns the next y force component
 * @return
 */
double circularForce::return_next_y() {
    return y_force(next_selector);
}

double circularForce::return_current_next_x() {
    return x_force(current_selector) + x_force(next_selector);
}

double circularForce::return_current_next_y() {
    return y_force(current_selector) + y_force(next_selector);
}

/**
 * @fn void circularForce::increment()
 * @brief increments the counter and checks if we have to change the selectors
 */
void circularForce::increment() {
    // look if we need to change the selectors
    selector_switchero();
    ++counter;
}

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
