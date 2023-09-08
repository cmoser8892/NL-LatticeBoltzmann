#ifndef NL_LATTICEBOLTZMANN_MARKER_H
#define NL_LATTICEBOLTZMANN_MARKER_H

#include "types.h"
#include "straight.h"
#include "functions.h"

/**
 * The marker class holds and distributes the markers on a surface for the IBM_OUTER method.
 * @note not sure yet how to handle finer grids prob fully dynamic ?!
 */
class markerPoints {
  private:
    double marker_distance = 0.75; /**< distance between markers on the surface */
    straightGenerator *sg = nullptr; /**< straights used for generation of the marker points */
    void walking(straight_t* s,double walker_dist);
  public:
    std::vector<point_t *> marker_points;
    explicit markerPoints(straightGenerator *s = nullptr, double md = 0.75);
    ~markerPoints();
    void distribute_markers();
    void distribute_markers_periodic(straight_t* line,boundaryType_t next, kernelType_t t);
    double return_marker_distance();
    void write_out_markers(bool write_file);
};

#endif // NL_LATTICEBOLTZMANN_MARKER_H
