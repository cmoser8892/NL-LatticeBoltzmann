#ifndef NL_LATTICEBOLTZMANN_MARKER_H
#define NL_LATTICEBOLTZMANN_MARKER_H

#include "types.h"
#include "straight.h"
#include "functions.h"

/**
 * The marker class holds and distributes the markers on a surface for the IBM_OUTER method.
 * @note not sure yet how to handle finer grids prob fully dynamic ?!
 */
class markerDistribution {
  private:
    double marker_distance = 0.75; /**< distance between markers on the surface */
    straightGenerator *sg = nullptr; /**< straights used for generation of the marker points */
  public:
    std::vector<point_t *> marker_points;
    explicit markerDistribution(straightGenerator *s = nullptr, double md = 0.75);
    ~markerDistribution();
    void distribute_markers();
    void individual_distribute_markers(straight_t s);
    double return_marker_distance();
    void write_out_markers(bool write_file);
    void set_marker_distance(double md);
};

#endif // NL_LATTICEBOLTZMANN_MARKER_H
