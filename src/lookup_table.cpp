
#include "lookup_table.h"
#include "helper_functions.h"
/// private
uint32_t Lookup::u_adc_converter(double u) {

}

void Lookup::fill_table() {

}

double Lookup::calculate_function(double cx, double cy, double ux, double uy) {

}

uint64_t Lookup::key_generation(uint32_t c_x, uint32_t c_y, uint32_t ux, uint32_t uy) {
    // bit interleaving of 1 bit, 1 bit, u_bit_rep, u_bit_rep

}

/// public
Lookup::Lookup(unsigned int how_many_u_bits,
               double floor,
               double ceiling,
               bool inter) :
u_bit_representation(how_many_u_bits),u_floor(floor),u_ceiling(ceiling),interpolation(inter){
    fill_table();
}

double Lookup::look_at_table(bool cx, bool cy, double ux, double uy) {
    double return_value = 0;
    uint32_t ux_bits = u_adc_converter(ux);
    uint32_t uy_bits = u_adc_converter(uy);
    uint64_t search_key = key_generation(cx,cy,ux_bits,uy_bits);
    // look in table
    if(auto found_iter = lookup_table.find(search_key); found_iter != lookup_table.end()) {
        return_value = found_iter->second;
    }
    // not in table recalculate
    else {
        non_hits++;
        return_value = calculate_function(double(cx),double(cy),ux,uy);
    }
    return return_value;
}