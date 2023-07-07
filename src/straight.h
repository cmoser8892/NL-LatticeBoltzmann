#ifndef NL_LATTICEBOLTZMANN_STRAIGHT_H
#define NL_LATTICEBOLTZMANN_STRAIGHT_H

#include "types.h"
#include "node.h"
#include "boundary_point_generator.h"
#include "helper_functions.h"
#define CRITICAL_RANGE 0 /**< The boundary intersection tests struggle when the mass center is near one */

/**
 * Straight definition.
 */
typedef struct straight {
    // s = p + t*d
    point_t point; /**<  Origin of the straight */
    vector_t direction; /**<  Direction of the straight line */
    // validity of the straight
    double max_t = 0; /**<  Length of that line */
}straight_t;

/**
 * Straight generator that connects the surface.
 */
class straightGenerator {
  private:
    // boundary points and pkh of the boundary points
    boundaryPointConstructor* points; /**< Pointer to the Boundary point */
    pointKeyHash full_pkh; /**<  Hash table of all the boundary points */
    std::vector<pointKeyHash*> pkhv; /**<  Hash table of each boundary structure */
    // global and individual mass centers
    point_t mass_center; /**<  mass_center of the boundary points to calculate intersection later */
    std::vector<point_t> individual_mass_centers; /**<  Individual mass centers of the boundary structures */
    // temporary object storages
    std::vector<straight_t *> temporary_creation; /**< temp 1 */
    std::vector<straight_t *> temporary_valid; /**<  temp 2 */
    std::vector<straight_t *> temporary; /**<  temp 3, yes we need 3 temporary storages to reform straight lines */
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
  public:
    // suface defined as middle point between two boundary points and a normal vector
    std::vector<straight_t *> surfaces; /**<  surfaces in a vector */
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
