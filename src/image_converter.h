
#ifndef NL_LATTICEBOLTZMANN_IMAGE_CONVERTER_H
#define NL_LATTICEBOLTZMANN_IMAGE_CONVERTER_H

#include "types.h"
#include "boundary_point_generator.h"
#include "helper_functions.h"
#include "bmp.h"

// definition of all the different parts of a bmp file
#define WHITE_COLOR_CODE_24_BIT     0xffffff     /**< The color of fluid nodes */
#define TARGET_SIZE_WINDOW          4            /**< Saved windows size */

/// converts a bmp image to a boundary
class imageConverter {
  private:
    std::unordered_map<colour_t ,colour_t> colors_used; /**< Map of the colors used in the image */
    pointKeyHash pkh; /**< Used to make the reformed points */
    int current_structure = 0; /**< Current structure counter */
    // functions
    void read();
    void detect_colors();
    void create_raw();
    void translate_reformed_into_structure();
    point_t update_position(point_t p);
  public:
    // public vars
    bmpReader* bmp_reader;
    rawPoints * raw = nullptr;
    boundaryPointConstructor* boundaries = nullptr;
    // functions
    explicit imageConverter(std::filesystem::path p);
    ~imageConverter();
    void init();
    void run();
    // individual tests of the run functionality
    void raw_run();
    void raw_cleanup();
    // helpers
    int return_number_of_colors();
    bool check_for_white_wet_nodes();
    unsigned long return_basic_size();
};

#endif // NL_LATTICEBOLTZMANN_IMAGE_CONVERTER_H
