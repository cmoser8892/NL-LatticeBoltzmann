//
// Created by christoph on 14.01.23.
//

#ifndef NL_LATTICEBOLTZMANN_STRAIGHT_H
#define NL_LATTICEBOLTZMANN_STRAIGHT_H

#include "types.h"
#include "boundary_point_generator.h"

typedef struct straight {
    point_t point;
    vector_t direction;
    double length = 0;
}straight_t;

class straight_generator {
  private:
    boundaryPointConstructor* points;
    point_t mass_center;
    void calculate_mass_center();
    void calculate_all_straights();
    void reduce();
  public:
    std::vector<straight_t*> straights;
    explicit straight_generator(boundaryPointConstructor* p);
    void init();
    void calculate_intersections();
};
#endif // NL_LATTICEBOLTZMANN_STRAIGHT_H
