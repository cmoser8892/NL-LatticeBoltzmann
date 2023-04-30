#ifndef NL_LATTICEBOLTZMANN_STRAIGHT_H
#define NL_LATTICEBOLTZMANN_STRAIGHT_H

#include "types.h"
#include "node.h"
#include "boundary_point_generator.h"
#include "helper_functions.h"

typedef struct straight {
    point_t point;
    vector_t direction;
}straight_t;

class straightGenerator {
  private:
    std::vector<pointKeyHash*> pkhv;
    boundaryPointConstructor* points;
    std::unordered_multimap<handle_t,handle_t> keys;
    point_t mass_center;
    std::vector<point_t> individual_mass_centers;
    void calculate_mass_center();
    void calculate_keys();
    void calculate_all_straights();
    int calculate_intersections(nodePoint_t* point);
  public:
    // suface defined as middle point between two boundary points and a normal vector
    std::vector<straight_t *> surfaces;
    explicit straightGenerator(boundaryPointConstructor* p);
    ~straightGenerator();
    void init();
    void init_test();
    bool node_inside(nodePoint_t *point);
    void delete_vector();
};
#endif // NL_LATTICEBOLTZMANN_STRAIGHT_H
