
#ifndef NL_LATTICEBOLTZMANN_IMAGE_CONVERTER_H
#define NL_LATTICEBOLTZMANN_IMAGE_CONVERTER_H

#include "types.h"
#include <filesystem>
#include <iostream>
#include <fstream>
#include "boundary_point_generator.h"

/// definition of all the different parts of a bmp file
// file header 14 byte
typedef struct __attribute__((packed)) BMPFileHeader {
    uint16_t file_type{0x4D42};     // stands for BM in hex
    uint32_t file_size{0};          // size in bytes
    uint16_t reserved_1{0};         // reserved always 0
    uint16_t reserved_2{0};         // reserved always 0
    uint32_t offset_data{0};        // start position of the pixel data
}BMPFileHeader_t;

// bitmap/info header 40 byte
typedef struct __attribute__((packed)) BMPInfoHeader {
    uint32_t size{0};                       // size of this header
    int32_t width{0};                       // width of the bitmap in pixel
    int32_t height{0};                      // height of the bitmap in pixel
    // depending on the sign bottom up or top down with different origins
    uint16_t planes{1};                     // always one
    uint16_t bit_count{0};                  // bits per pixel
    uint32_t compression{0};                // 0 (24b) or 3 (32b) only uncompressed images
    uint32_t size_image{0};                 // 0 for uncompressed images
    int32_t x_pixels_per_meter{0};          // horizontal resolution of the image
    int32_t y_pixels_per_meter{0};          // vertical resolution of the image
    uint32_t colors_used{0};                // numbers of colours used
    uint32_t colors_important{0};           // if 0 all the colors are important otherwise that number
}BMPInfoHeader_t;

// color table
typedef struct __attribute__((packed)) BMPColorTable {
    uint32_t red_mask{0x00ff0000};         // Bit mask for the red channel
    uint32_t green_mask{0x0000ff00};       // Bit mask for the green channel
    uint32_t blue_mask{0x000000ff};        // Bit mask for the blue channel
    uint32_t alpha_mask{0xff000000};       // Bit mask for the alpha channel
    uint32_t color_space_type{0x73524742}; // Default "sRGB" (0x73524742)
    uint32_t unused[16]{ 0 };              // Unused data for sRGB color space
}BMPColorTable_t;

// definition of all the data fields
typedef struct BMP {
    BMPFileHeader_t file_header;
    BMPInfoHeader_t info_header;
    BMPColorTable_t color_header;
    std::vector<uint8_t> data;
}BMP_t;

/// converts a bmp image to a boundary
class imageConverter {
  private:
    std::filesystem::path path; // path to the the bmp file
    boundaryPointConstructor* boundaries;
    std::unordered_map<colour_t, boundaryType_t> mapping;
    BMP_t bmp;
    // functions
    void read();
    uint32_t make_stride_aligned(uint32_t align_stride,uint32_t row_stride);
  public:
    imageConverter(std::filesystem::path p);
    void init();
    void map_colours_to_boundaries();
    void communicate_colour_decision();
};

#endif // NL_LATTICEBOLTZMANN_IMAGE_CONVERTER_H
