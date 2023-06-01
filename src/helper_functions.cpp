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
    if(out.is_open()) {
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
    else {
        throw std::runtime_error("Open file has failed");
    }
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
 * @fn double calculate_distance(point_t* p1, point_t* p2)
 * @brief calculate the distance between two points
 * @param p1
 * @param p2
 * @return
 */
double calculate_distance(point_t* p1, point_t* p2) {
    // hide this ugliness
    return (*p1 - *p2).norm();
}

/**
 * @fn int conical_delta(int a, int b)
 * @brief returns 1 if a == b otherwise 0
 * @param a
 * @param b
 * @return
 */
int conical_delta(int a, int b) {
    return a == b;
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
