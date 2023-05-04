#ifndef NL_LATTICEBOLTZMANN_STRAIGHT_H
#define NL_LATTICEBOLTZMANN_STRAIGHT_H

#include "types.h"
#include "node.h"
#include "boundary_point_generator.h"
#include "helper_functions.h"

typedef struct straight {
    // s = p + t*d
    point_t point;
    vector_t direction;
    // validity of the straight
    double min_t = 0;
    double max_t = 0;
}straight_t;

class straightGenerator {
    /**
     * note: surfaces are the boundary point to the next boundary point
     * connected via a direction vector
     * it is a straight line
     */
  private:
    // boundary points and pkh of the boundary points
    boundaryPointConstructor* points;
    std::vector<pointKeyHash*> pkhv;
    // global and individual mass centers
    point_t mass_center;
    std::vector<point_t> individual_mass_centers;
    // temporary object storages
    std::vector<straight_t *> temporary_creation;
    std::vector<straight_t *> temporary_valid;
    // functions
    void calculate_mass_center();
    void calculate_keys();
    void calculate_all_straights();
    int calculate_intersections(nodePoint_t* point);
    void straight_create(int bs);
    void straight_self_test(int bs);
    bool straight_closer_test(int bs, straight_t* self, straight_t* partner);
    bool straight_better_candidate_test(straight_t* candidate, straight* partner);
    void straight_set_t_values(int bs);
    double go_through_vector(int bs, straight_t* self,int direction);
    void straight_test_creation(int bs);
  public:
    // suface defined as middle point between two boundary points and a normal vector
    std::vector<straight_t *> surfaces;
    explicit straightGenerator(boundaryPointConstructor* p);
    ~straightGenerator();
    void init();
    void init_test();
    bool node_inside(nodePoint_t *point);
    void delete_vector();
    void delete_keys();
};

double calculate_intersection(straight_t* ray, straight_t* surface);
#endif // NL_LATTICEBOLTZMANN_STRAIGHT_H
