
#ifndef NL_LATTICEBOLTZMANN_BOUNDARY_POINT_GENERATOR_H
#define NL_LATTICEBOLTZMANN_BOUNDARY_POINT_GENERATOR_H

#include "types.h"

/**
 * class to generate/hold boundary points from other structures
 * first ingredient to the init
 */

typedef struct boundaryPoint {
    point_t point;
    boundaryType_t type;
}boundaryPoint_t;

class boundaryPointConstructor {
  private:
    // nada
  public:
    // holds all the boundary points
    std::vector<boundaryPoint_t *> boundary_points;
    // sizes and limits of the construction field
    point_t size;
    point_t limits;
    //
    explicit boundaryPointConstructor(point_t s);
    // initizlies a quader
    void one_direction(int limit,vector_t dir,point_t* start, boundaryType_t b);
    void init_quader();
};

#endif // NL_LATTICEBOLTZMANN_BOUNDARY_POINT_GENERATOR_H
