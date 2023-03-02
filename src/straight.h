//
// Created by christoph on 14.01.23.
//

#ifndef NL_LATTICEBOLTZMANN_STRAIGHT_H
#define NL_LATTICEBOLTZMANN_STRAIGHT_H

#include "types.h"
#include "node.h"
#include "boundary_point_generator.h"

typedef struct straight {
    point_t point;
    vector_t direction;
}straight_t;

// for easier nomenclature
using surface_t = straight_t;

class straightGenerator {
  private:
    boundaryPointConstructor* points;
    point_t mass_center;
    void calculate_mass_center();
    void calculate_all_straights();
    int calculate_intersections(nodePoint_t* point);
  public:
    // suface defined as middle point between two boundary points and a normal vector
    std::vector<surface_t *> surfaces;
    explicit straightGenerator(boundaryPointConstructor* p);
    ~straightGenerator();
    void init();
    bool node_inside(nodePoint_t *point);
    void delete_vector();
};
#endif // NL_LATTICEBOLTZMANN_STRAIGHT_H
