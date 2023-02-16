//
// Created by christoph on 16.02.23.
//

#ifndef NL_LATTICEBOLTZMANN_LOOKUP_TABLE_H
#define NL_LATTICEBOLTZMANN_LOOKUP_TABLE_H

#include "types.h"

class Lookup {
  private:
    // how many bits are used for the representation of ux and uy
    uint32_t u_bit_representation = 0;
    double u_floor;
    double u_ceiling;
    bool interpolation;
    std::unordered_map<unsigned, double> lookup_table;
    // convert a double to the representation
    uint32_t u_adc_converter(double u);
    void fill_table();
    double calculate_function(double cx, double cy, double ux, double uy);
    uint64_t key_generation(uint32_t c_x, uint32_t c_y, uint32_t ux, uint32_t uy);
  public:
    uint64_t non_hits = 0;
    Lookup(unsigned how_many_u_bits, double u_floor, double u_ceiling, bool interpolation);
    double look_at_table(bool cx, bool cy, double ux, double uy);
};

#endif // NL_LATTICEBOLTZMANN_LOOKUP_TABLE_H
