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
    double max_t = 0;
}straight_t;

class straightGenerator {
    /**
     * note: surfaces are the boundary point to the next boundary point
     * connected via a direction vector
     * it is a straight line
     */
  private:
    int intersection_test_alternate_state_1 = 0;
    int intersection_test_alternate_state_2 = 0;
    // boundary points and pkh of the boundary points
    boundaryPointConstructor* points;
    pointKeyHash full_pkh;
    std::vector<pointKeyHash*> pkhv;
    // global and individual mass centers
    point_t mass_center;
    std::vector<point_t> individual_mass_centers;
    // temporary object storages
    std::vector<straight_t *> temporary_creation;
    std::vector<straight_t *> temporary_valid;
    std::vector<straight_t *> temporary;
    // functions
    void calculate_mass_center();
    void detect_boundary_proximity_main_mass_center();
    void calculate_keys();
    void calculate_all_straights(); // old simpler method only works for concave surfaces with no bumps or anything
    int calculate_intersections(point_t node_point, point_t *individual_mc);
    int calculate_intersections_redundant(nodePoint_t* point);
    bool calculate_intersections_star_node_point(nodePoint_t* point);
    // creation related methods
    void straight_create(int bs);
    void straight_reduce(int bs);
    void find_surface_boundary_points(int bs);
    void look_for_bumps(int bs);
    void temporary_to_surface(int bs);
    void post_run_messages();
  public:
    // suface defined as middle point between two boundary points and a normal vector
    std::vector<straight_t *> surfaces;
    explicit straightGenerator(boundaryPointConstructor* p);
    ~straightGenerator();
    void init();
    bool node_inside_simple(nodePoint_t *point);
    bool node_inside_star(nodePoint_t *point);
    void delete_vector();
    void delete_keys();
    void write_out_surface();
};

bool straight_better_candidate_test(straight_t* candidate, straight* partner);
double calculate_intersection(straight_t* ray, straight_t* surface);
bool compare_bumps_sort(const std::tuple<double,double,handle_t,straight_t*> &a,
                        const std::tuple<double,double,handle_t,straight_t*> &b);

// todo fix errors related to mass center placement in the check inside algorithm :)
// also todo  make 3 mcs for a more stable algorithm
// alternatively add a surface identifier to a surface if we cut the same surface more than once we know that we got sth convex
// add more mass centers based on
// how do i test sth to be convex?!
// write one more test related to the construction oposing bumps should not be a prob thou
// min out + general cleanup
#endif // NL_LATTICEBOLTZMANN_STRAIGHT_H
