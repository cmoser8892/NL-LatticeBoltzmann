#ifndef NL_LATTICEBOLTZMANN_BOUNDARY_POINT_GENERATOR_H
#define NL_LATTICEBOLTZMANN_BOUNDARY_POINT_GENERATOR_H

#include "types.h"
#include "helper_functions.h"
/**
 * class to generate/hold boundary points from other structures
 * first ingredient to the init
 */

// definition of a general boundary point
typedef struct boundaryPoint {
    handle_t h; // the handle is forgotten when moved to node_infos..., it gets a new on
    point_t point;
    nodeIdentifier_t dw;
    boundaryType_t type;
}boundaryPoint_t;

// raw unordered point clouds
class rawBoundaryPoints {
    // enum return codes for boarder-checker
    // only visible in the class functions no need otherwise
    typedef enum border_return_code {
        INSIDE = 0,
        BOARDER,
        CORNER
    }border_return_code_t;
  private:
    // reminder that valid handles start at 1
    handle_t current_handle = 0;
    point_t size;
    point_t limits;
    pointKeyHash pkh;
    void fill_keys();
    border_return_code_t check_boarder(boundaryPoint_t &b);
    int set_max_neighbors(border_return_code_t b);
  public:
    std::vector<boundaryPoint_t*> raw_boundary_points;
    std::vector<boundaryPoint_t*> reformed_boundary_points;
    rawBoundaryPoints(point_t size);
    ~rawBoundaryPoints();
    void delete_raw_boundary_points();
    void delete_reformed_boundary_points();
    void rewrite_reformed_boundary_handles();
    void visualize_2D_boundary();
    void read_in_bounce_back(point_t p);
    void read_in_bounce_back(coordinate_t p);
    void reduce();
};

// Description of one closed surface/structure
class boundaryStructure {
  public:
    ~boundaryStructure();
    void rewrite_reformed_boundary_handles();
    // data
    std::vector<boundaryPoint_t*> boundary_points;
};

class boundaryPointConstructor {
  private:
    // reminder valid handles start at 1
    handle_t added_handle = 0;

  public:
    // holds all the boundary points
    long current_structure = -1;
    std::vector<boundaryStructure *> boundary_structures;
    // sizes and limits of the construction field
    point_t size;
    point_t limits;
    //
    explicit boundaryPointConstructor(point_t s);
    ~boundaryPointConstructor();
    // initializes a general structure
    void init_structure();
    // init points and rewrite handles
    void set_point(point_t *p, boundaryType_t b);
    void set_point(handle_t h, point_t *p, boundaryType_t b);
    void rewrite_handles();
    void delete_structures();
    // visualize the resulting structure
    void visualize_2D_boundary();
    // total nodes
    long total_boundary_nodes();
    // different functions to set up different simple structures (heavily used in testing)
    void one_direction(int limit, vector_t dir, point_t *start, boundaryType_t b);
    void steps_direction(int steps, vector_t dir, point_t *start, boundaryType_t b);
    void corner_creation(vector_t dir, point_t *start, boundaryType_t b);
    void init_quader();
    void init_chopped_quader(point_t start, int devider);
    void init_quader(point_t start);
    void init_quader(point_t start, vector_t size);
    void init_sliding_lid();
    void init_chopped_sliding_lid(point_t start, int devider);
    void init_quader_side_chopped(point_t start, int chopsize);
    void init_sliding_lid_side_chopped(point_t start, int chopsize);
    void init_sliding_lid_inner(point_t start, point_t cntiues, vector_t inner_size);
    void init_poiseuille_flow();
    void pressure_inlet();
};

bool sorter_wet_dry_boundaries(boundaryPoint_t * p1, boundaryPoint_t * p2);
#endif // NL_LATTICEBOLTZMANN_BOUNDARY_POINT_GENERATOR_H
