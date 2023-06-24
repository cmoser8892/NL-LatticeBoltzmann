#ifndef NL_LATTICEBOLTZMANN_BOUNDARY_POINT_GENERATOR_H
#define NL_LATTICEBOLTZMANN_BOUNDARY_POINT_GENERATOR_H

#include "types.h"
#include "helper_functions.h"
#include "helper_classes.h"

/**
 * Contains the information of a boundary point (handle, point, dry/wet, boundary type).
 */
typedef struct boundaryPoint {
    handle_t h; /**< handle to the object, only unique in the boundary point cloud, will get a new one in nodes*/
    point_t point; /**< position of the boundary point */
    nodeIdentifier_t dw; /**< dry or wet boundary */
    boundaryType_t type; /**< boundary type (bounce back, periodic, etc.*/
}boundaryPoint_t;

/**
 * Holds all possible points that are not part of a fluid, the class reduce those points just to boundary points on the boarder
 */
class rawPoints {
    /**
     * Description of cases where a point is, inside on the boarder or at a corner.
     */
    typedef enum border_return_code {
        INSIDE = 0, /**< inside of the bulk of the raw points, means I have 8 neighbors */
        BOARDER, /**< at the boarder of the bulk of raw points */
        CORNER /**< special boarder case of a boundary */
    }border_return_code_t;
  private:
    handle_t current_handle = 0; /**< counter of the currently used handles to be given */
    point_t size; /**< overall size of the object */
    point_t limits; /**< size -1 */
    pointKeyHash pkh; /**< hashing class used to find the immediate neighbors of a point */
    void fill_keys();
    border_return_code_t check_boarder(boundaryPoint_t &b);
    int set_max_neighbors(border_return_code_t b);
    int set_min_neighbors(border_return_code_t b);
    bool judge_add_up_found_velocities_vector(vector_t a);
  public:
    std::vector<boundaryPoint_t*> raw_boundary_points; /**< vector of all boundary points */
    std::vector<boundaryPoint_t*> reformed_boundary_points; /**< vector of the variable boundary points */
    explicit rawPoints(point_t size);
    ~rawPoints();
    void delete_raw_boundary_points();
    void delete_reformed_boundary_points();
    void rewrite_reformed_boundary_handles();
    void visualize_2D_boundary();
    void read_in_bounce_back(point_t p);
    void read_in_bounce_back(coordinate_t p);
    void reduce();
};

/**
 * This class holds a closed amount of boundary points to be used when constructing surfaces.
 */
class boundaryStructure {
  public:
    ~boundaryStructure();
    void rewrite_reformed_boundary_handles();
    std::vector<boundaryPoint_t*> boundary_points; /**< data vector holder of boundary points */
};
/**
 * This class can either construct a boundary with a set of functions or process reduced raw points into a structured boundary point cloud.
 */
class boundaryPointConstructor {
  private:
    handle_t added_handle = 0;
  public:
    // holds all the boundary points
    long current_structure = -1;
    std::vector<boundaryStructure *> boundary_structures;
    // sizes and limits of the construction field
    point_t size;
    point_t limits;
    // functions
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
    void init_chopped_quader(point_t start,point_t size, int devider);
    void init_quader(point_t start, vector_t size);
    void init_sliding_lid();
    void init_chopped_sliding_lid(point_t start, point_t size ,int devider);
    void init_quader_side_chopped(point_t start, int chopsize);
    void init_sliding_lid_side_chopped(point_t start, int chopsize);
    void init_sliding_lid_inner(point_t start,vector_t outer_size, point_t cntiues, vector_t inner_size);
    void init_poiseuille_flow();
    void pressure_inlet();
};

bool sorter_wet_dry_boundaries(boundaryPoint_t * p1, boundaryPoint_t * p2);
#endif // NL_LATTICEBOLTZMANN_BOUNDARY_POINT_GENERATOR_H
