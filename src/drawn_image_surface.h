#ifndef NL_LATTICEBOLTZMANN_DRAWN_IMAGE_SURFACE_H
#define NL_LATTICEBOLTZMANN_DRAWN_IMAGE_SURFACE_H

#include "bmp.h"
#include "types.h"
#include "functions.h"
#include "helper_classes.h"
#include "straight.h"

class surfaceDrawer {
  private:
    std::filesystem::path path;
    void add_surface(point_t current, point_t next);
  public:
    straightGenerator surface_storage;
    // functions
    explicit surfaceDrawer(std::filesystem::path p);
    ~surfaceDrawer();
    void run();
};

#endif // NL_LATTICEBOLTZMANN_DRAWN_IMAGE_SURFACE_H
