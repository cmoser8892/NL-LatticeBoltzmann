#ifndef NL_LATTICEBOLTZMANN_MARKER_H
#define NL_LATTICEBOLTZMANN_MARKER_H

#include "types.h"
#include "straight.h"
#include "functions.h"
#include "nodeGenerator.h"

/**
 * The marker class holds and distributes the markers on a surface for the IBM method.
 * @note not sure yet how to handle finer grids prob fully dynamic ?!
 * @todo not yet sure on the dual role of the markers
 */
class markerIBM {
  private:
    int kernel_length = 2;         /**< search space of the kernel */
    double marker_distance = 0.75; /**< distance between markers on the surface */
    straightGenerator *sg = nullptr;

  public:
    markerIBM(straightGenerator *s = nullptr, int kd = 2, double md = 0.75);
    std::vector<point_t *> marker_points;
    void find_neighborhood();
    void distribute_markers();
};

#endif // NL_LATTICEBOLTZMANN_MARKER_H
