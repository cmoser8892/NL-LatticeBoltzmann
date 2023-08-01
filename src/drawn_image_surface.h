#ifndef NL_LATTICEBOLTZMANN_DRAWN_IMAGE_SURFACE_H
#define NL_LATTICEBOLTZMANN_DRAWN_IMAGE_SURFACE_H

#include "bmp.h"
#include "types.h"
#include "functions.h"

class surfaceDrawer {
  private:
    void read();
    void convert_points();
  public:
    // public vars
    bmpReader* bmp_reader;
    std::vector<point_t> points;
    // functions
    explicit surfaceDrawer(std::filesystem::path p);
    ~surfaceDrawer();
    void init();
    void run();
};

#endif // NL_LATTICEBOLTZMANN_DRAWN_IMAGE_SURFACE_H
