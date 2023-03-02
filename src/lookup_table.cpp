#include "lookup_table.h"
#include "helper_functions.h"
/// private
void lookup::check_u_value_outside(double u) {
    if(u < u_floor) {
        error_flag = true;
    }
    if(u> u_ceiling) {
        error_flag = true;
    }
}

double lookup::u_adc_number(double u) {
    double ret = (u-u_floor)/u_step;
    return ret;
}

uint32_t lookup::u_adc_converter_lower(double u) {
    return uint32_t(std::floor(u_adc_number(u)));
}

uint32_t lookup::u_adc_converter_higher(double u) {
    return uint32_t(std::ceil(u_adc_number(u)));
}

void lookup::fill_table() {
    double ux = u_floor;
    double uy = u_floor;
    // enumerate the values and generate bit representations
    for(int i = -1; i < 2; ++i) {
        auto cx_bit = reduce_32_2(i);
        for(int j = -1; j < 2; j++) {
            auto cy_bit = reduce_32_2(j);
            for(int k = 0; k < u_bit_number; ++k) {
                ux = u_floor +k*u_step;
                uint32_t ux_bit = u_adc_converter_lower(ux);
                for(int l = 0; l < u_bit_number; ++l) {
                    uy = u_floor + l* u_step;
                    uint32_t uy_bit = u_adc_converter_lower(uy);
                    // calculate and put key + number together
                    uint64_t key = key_generation(cx_bit,cy_bit,ux_bit,uy_bit);
                    double value = calculate_function(i,j,ux,uy);
                    lookup_table.emplace(key, value);
                }
            }
        }
    }
}

double lookup::calculate_function(double cx, double cy, double ux, double uy) {
    return (1 + 3*cx*ux + 3*cy*uy + 4.5*cx*cx*ux*ux + 9*cx*ux*cy*uy + 4.5*cy*cy*uy*uy - 1.5*ux*ux - 1.5*uy*uy);
}

uint64_t lookup::key_generation(uint32_t cx, uint32_t cy, uint32_t ux, uint32_t uy) {
    // bit interleaving of 1 bit, 1 bit, u_bit_rep, u_bit_rep
    // key looks ux_bit interleaved with uy_bit cy_bit cx_bit
    // masking cx and cy (should be 0x1 or 0x0 anyways)
    uint64_t return_value = bit_interleaving_2d(ux,uy);
    // todo possible to reduce the bit code by one got the truthtables how do i reduce them again
    return_value |= (((cx << 2)| cy) << u_bit_representation) << u_bit_representation;
    return return_value;
}

/// public
lookup::lookup(unsigned int how_many_u_bits,
               double floor,
               double ceiling,
               bool inter) :
u_bit_representation(how_many_u_bits),u_floor(floor),u_ceiling(ceiling),interpolation(inter){
    u_bit_number = 1 << how_many_u_bits;
    u_step = (u_ceiling-u_floor)/u_bit_number;
    fill_table();
}

void lookup::set_bypass(bool b) {
    bypass = b;
}

double lookup::look_at_table(int cx, int cy, double ux, double uy) {
    double return_value;
    if(bypass) {
        return calculate_function(cx,cy,ux,uy);
    }
    // check if were even in the bracket
    check_u_value_outside(ux);
    check_u_value_outside(uy);
    // check for errors in the conversion to pure bits key generation otherwise wont work
    if(!error_flag) {
        uint32_t ux_bits = u_adc_converter_lower(ux);
        uint32_t uy_bits = u_adc_converter_lower(uy);
        uint32_t cx_bits = reduce_32_2(cx);
        uint32_t cy_bits = reduce_32_2(cy);
        uint64_t search_key = key_generation(cx_bits, cy_bits, ux_bits, uy_bits);
        // look in table
        if (auto found_iter = lookup_table.find(search_key); found_iter != lookup_table.end()) {
            ++table_hits;
            return_value = found_iter->second;
        }
        // not in table recalculate
        else {
            ++non_hits;
            return_value = calculate_function(double(cx), double(cy), ux, uy);
        }
    }
    else {
        ++recalc;
        error_flag = false;
        return_value = calculate_function(double(cx),double(cy),ux,uy);
    }
    return return_value;
}

double lookup::look_at_table(uint32_t cx_bit, uint32_t cy_bit, uint32_t ux_bit, uint32_t uy_bit) {
    double return_value = 0;
    uint64_t search_key = key_generation(cx_bit, cy_bit, ux_bit, uy_bit);
    // look in table
    if (auto found_iter = lookup_table.find(search_key); found_iter != lookup_table.end()) {
        ++table_hits;
        return_value = found_iter->second;
    }
    // not in table recalculate
    else {
        ++non_hits;

    }
    return return_value;
}

array_t lookup::equilibrium(node* n) {
    array_t return_array;
    return_array.setZero(CHANNELS);
    double ux = n->u(0);
    double uy = n->u(1);
    double rho = n->rho;
    bool recalculate = false;
    if(!bypass){
        check_u_value_outside(ux);
        check_u_value_outside(uy);
        // calculate weather or not inside
        if (!error_flag) {
            uint32_t ux_bit = u_adc_converter_lower(ux);
            uint32_t uy_bit = u_adc_converter_lower(uy);
            for (int i = 0; i < CHANNELS; ++i) {
                uint32_t cx_bit = reduce_32_2(velocity_set.col(i).x());
                uint32_t cy_bit = reduce_32_2(velocity_set.col(i).y());
                double w = weights.col(i).x();
                return_array(i) = w * rho * look_at_table(cx_bit, cy_bit, ux_bit, uy_bit);
            }
        }
        else {
            recalculate = true;
        }
    }
    else {
        recalculate = true;
    }
    if(recalculate) {
        error_flag = false;
        for (int i = 0; i < CHANNELS; ++i) {
            ++recalc;
            double cx = velocity_set.col(i).x();
            double cy = velocity_set.col(i).y();
            double w = weights.col(i).x();
            return_array(i) = w * rho * calculate_function(cx, cy, ux, uy);
        }
    }
    return return_array;
}