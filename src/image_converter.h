
#ifndef NL_LATTICEBOLTZMANN_IMAGE_CONVERTER_H
#define NL_LATTICEBOLTZMANN_IMAGE_CONVERTER_H

#include "types.h"
#include <filesystem>
#include "boundary_point_generator.h"

/// converts a bmp image to a boundary structure
class imageConverter {
  private:
    std::filesystem::path path;
    boundaryPointConstructor* boundaries;
    std::unordered_map<colour_t, boundaryType_t> mapping;
  public:
    imageConverter(std::filesystem::path p);
    void check_correct_format();
    void map_colours_to_boundaries();
    void communicate_colour_decision();
};


#endif // NL_LATTICEBOLTZMANN_IMAGE_CONVERTER_H
