
#include "lookup_table.h"
#include "helper_functions.h"
/// private
double Lookup::u_adc_number(double u) {
    double ret = (u-u_floor)*u_bit_number;
    // checks if within bounds
    if(ret < u_floor) {
        error_flag = true;
    }
    if(ret > u_floor) {
        error_flag = true;
    }
    return ret;
}
uint32_t Lookup::u_adc_converter_lower(double u) {
    return uint32_t(std::floor(u_adc_number(u)));
}

uint32_t Lookup::u_adc_converter_higher(double u) {
    return uint32_t(std::ceil(u_adc_number(u)));
}

void Lookup::fill_table() {
    double cx, cy = 0;
    double ux = u_floor;
    double uy = u_floor;
    // enumerate the values and generate bit representations
    for(int i; i < 2; ++i) {
        auto cx_bit = uint32_t(cx);
        for(int j; j < 2; j++) {
            auto cy_bit = uint32_t(cy);
            for(int k; k < u_bit_number; ++k) {
                uint32_t ux_bit = u_adc_converter_lower(ux);
                for(int l; l < u_bit_number; ++l) {
                    uint32_t uy_bit = u_adc_converter_lower(uy);
                    // calculate and put key + number together
                    uint64_t key = key_generation(cx_bit,cy_bit,ux_bit,uy_bit);
                    double value = calculate_function(cx,cy,ux,uy);
                    lookup_table.emplace(key, value);
                    // functional ribcage
                    uy += u_step;
                }
                ux += u_step;
            }
            cy+=1;
        }
        cx+=1;
    }
}

double Lookup::calculate_function(double cx, double cy, double ux, double uy) {
    return 2 + 6*cx*ux + 6*cy*uy + 9*cx*cx*ux*ux + 18*cx*ux*cy*uy + 9*cy*cy*uy*uy - 3*ux*ux -3*uy*uy;
}

uint64_t Lookup::key_generation(uint32_t cx, uint32_t cy, uint32_t ux, uint32_t uy) {
    // bit interleaving of 1 bit, 1 bit, u_bit_rep, u_bit_rep
    // key looks ux_bit interleaved with uy_bit cy_bit cx_bit
    // masking cx and cy (should be 0x1 or 0x0 anyways)
    uint64_t return_value = cx & 0x1;
    cy &= 0x1;
    return_value &= cy << 1;
    return_value &= bit_interleaving_2d(ux,uy) << 2;
    return return_value;
}

/// public
Lookup::Lookup(unsigned int how_many_u_bits,
               double floor,
               double ceiling,
               bool inter) :
u_bit_representation(how_many_u_bits),u_floor(floor),u_ceiling(ceiling),interpolation(inter){
    u_bit_number = 1 << how_many_u_bits;
    u_step = (u_ceiling-u_floor)/u_bit_number;
    fill_table();
}

void Lookup::set_bypass(bool b) {
    bypass = b;
}

double Lookup::look_at_table(bool cx, bool cy, double ux, double uy) {
    double return_value = 0;
    if(bypass) {
        return calculate_function(cx,cy,ux,uy);
    }
    uint32_t ux_bits = u_adc_converter_lower(ux);
    uint32_t uy_bits = u_adc_converter_lower(uy);
    // check for errors in the conversion to pure bits key generation otherwise wont work
    if(!error_flag) {
        uint64_t search_key = key_generation(cx, cy, ux_bits, uy_bits);
        // look in table
        if (auto found_iter = lookup_table.find(search_key); found_iter != lookup_table.end()) {
            return_value = found_iter->second;
        }
        // not in table recalculate
        else {
            ++non_hits;
            return_value = calculate_function(double(cx), double(cy), ux, uy);
        }
    }
    else {
        error_flag = false;
        return_value = calculate_function(double(cx),double(cy),ux,uy);
    }
    return return_value;
}