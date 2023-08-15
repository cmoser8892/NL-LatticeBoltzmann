#ifndef NL_LATTICEBOLTZMANN_DRAWN_IMAGE_SURFACE_H
#define NL_LATTICEBOLTZMANN_DRAWN_IMAGE_SURFACE_H

#include "bmp.h"
#include "types.h"
#include "functions.h"
#include "helper_classes.h"
#include "straight.h"

/**
 * The surface Drawer class hides all the opencv stuff
 */
class surfaceDrawer {
  private:
    std::filesystem::path path;  /**< Path to the image */
    void add_surface(point_t current, point_t next);
    bool selector(std::vector<int> sel, int current_contour);
  public:
    straightGenerator surface_storage; /**< The surface generator we write our surface to */
    // functions
    explicit surfaceDrawer(std::filesystem::path p);
    void run();
    void run_selective(std::vector<int> selector);
    void run_non_connecting(std::vector<int> selector, bool b);
    void close_open_surface(vector_t draw_size);
};

#endif // NL_LATTICEBOLTZMANN_DRAWN_IMAGE_SURFACE_H
