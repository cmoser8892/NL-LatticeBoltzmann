#ifndef NL_LATTICEBOLTZMANN_DRAWN_IMAGE_SURFACE_H
#define NL_LATTICEBOLTZMANN_DRAWN_IMAGE_SURFACE_H

#include "bmp.h"
#include "types.h"
#include "functions.h"
#include "src_foreign/spline.h"

class surfaceDrawer {
  private:
    void read();
    void convert_points();
    void spline_application();
  public:
    // public vars
    bmpReader* bmp_reader;
    // raw points
    std::vector<point_t> points;
    //
    // functions
    explicit surfaceDrawer(std::filesystem::path p);
    ~surfaceDrawer();
    void init();
    void run();
};

#endif // NL_LATTICEBOLTZMANN_DRAWN_IMAGE_SURFACE_H
