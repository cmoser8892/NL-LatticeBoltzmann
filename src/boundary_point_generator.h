#ifndef NL_LATTICEBOLTZMANN_BOUNDARY_POINT_GENERATOR_H
#define NL_LATTICEBOLTZMANN_BOUNDARY_POINT_GENERATOR_H

#include "types.h"

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

// Description of one closed surface/structure
class boundaryStructure {
  public:
    ~boundaryStructure();
    //
    std::vector<boundaryPoint_t *> boundary_points;

};

class boundaryPointConstructor {
  private:
    // nada
    handle_t added_handle = 0;
  public:
    // holds all the boundary points
    long current_structure = -1;
    std::vector<boundaryStructure*> boundary_structures;
    // sizes and limits of the construction field
    point_t size;
    point_t limits;
    //
    explicit boundaryPointConstructor(point_t s);
    ~boundaryPointConstructor();
    // initizlies a quader
    void init_structure();

    void one_direction(int limit,vector_t dir,point_t* start, boundaryType_t b);
    void steps_direction(int steps, vector_t dir, point_t* start, boundaryType_t b);
    void corner_creation(vector_t dir, point_t* start, boundaryType_t b);
    void set_point(point_t* p, boundaryType_t b);
    void init_quader();
    void init_chopped_quader(point_t start, int devider);
    void init_quader(point_t start);
    void init_quader(point_t start,vector_t size);
    void init_sliding_lid();
    void init_chopped_sliding_lid(point_t start, int devider);
    void init_quader_side_chopped(point_t start, int chopsize);
    void init_sliding_lid_side_chopped(point_t start, int chopsize);
    void init_sliding_lid_inner(point_t start, point_t cntiues, vector_t inner_size);
    void init_poiseuille_flow();
    void pressure_inlet();
    void delete_structures();
    // visualize
    void visualize_2D_boundary(int size);
    // total nodes
    long total_boundary_nodes();
};

bool sorter_wet_dry_boundaries(boundaryPoint_t * p1, boundaryPoint_t * p2);
#endif // NL_LATTICEBOLTZMANN_BOUNDARY_POINT_GENERATOR_H
