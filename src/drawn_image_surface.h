#ifndef NL_LATTICEBOLTZMANN_DRAWN_IMAGE_SURFACE_H
#define NL_LATTICEBOLTZMANN_DRAWN_IMAGE_SURFACE_H

#include "bmp.h"
#include "types.h"
#include "functions.h"
#include "helper_classes.h"
#include "straight.h"

class surfaceDrawer {
  private:
    handle_t handle_runner = 0;
    double step = 0.25;
    double range = 2;
    rangingPointKeyHash rpkh;
    void read();
    void convert_points();
    void fill_hashtable();
    point_t interpolate_around(point_t p);
    vector_t determine_init_surface_direction(point_t p);
    bool look_for_last(point_t current);
    void add_surface(point_t current, point_t previous);
  public:
    // public vars
    bmpReader* bmp_reader;
    // raw points
    std::vector<point_t> points;
    //
    straightGenerator surface_storage;
    // functions
    explicit surfaceDrawer(std::filesystem::path p);
    ~surfaceDrawer();
    void init();
    void run();
};

#endif // NL_LATTICEBOLTZMANN_DRAWN_IMAGE_SURFACE_H
