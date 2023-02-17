//
// Created by christoph on 16.02.23.
//

#ifndef NL_LATTICEBOLTZMANN_LOOKUP_TABLE_H
#define NL_LATTICEBOLTZMANN_LOOKUP_TABLE_H

#include "types.h"
#include "node.h"

class lookup {
  private:
    // how many bits are used for the representation of ux and uy
    uint32_t u_bit_representation = 0;
    uint32_t u_bit_number;
    double u_floor;
    double u_ceiling;
    double u_step;
    bool interpolation;
    bool error_flag = false;
    bool bypass = false;
    std::unordered_map<unsigned, double> lookup_table;
    // convert a double to the representation
    double u_adc_number(double u);
    uint32_t u_adc_converter_lower(double u);
    uint32_t u_adc_converter_higher(double u);
    void fill_table();
    double calculate_function(double cx, double cy, double ux, double uy);
    uint64_t key_generation(uint32_t c_x, uint32_t c_y, uint32_t ux, uint32_t uy);
  public:
    uint64_t non_hits = 0;
    uint64_t table_hits = 0;
    uint64_t recalc = 0;
    lookup(unsigned how_many_u_bits, double u_floor, double u_ceiling, bool interpolation);
    void set_bypass(bool);
    double look_at_table(int cx, int cy, double ux, double uy);
    array_t equilibrium(node* n);
};

#endif // NL_LATTICEBOLTZMANN_LOOKUP_TABLE_H
