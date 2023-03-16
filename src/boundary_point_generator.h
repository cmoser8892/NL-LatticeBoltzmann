#ifndef NL_LATTICEBOLTZMANN_BOUNDARY_POINT_GENERATOR_H
#define NL_LATTICEBOLTZMANN_BOUNDARY_POINT_GENERATOR_H

#include "types.h"

/**
 * class to generate/hold boundary points from other structures
 * first ingredient to the init
 */

// definition of a general boundary point
typedef struct boundaryPoint {
    point_t point;
    boundaryType_t type;
}boundaryPoint_t;

// Description of one closed surface/structure
class boundaryStructure {
  public:
    ~boundaryStructure();
    std::vector<boundaryPoint_t *> boundary_points;

};

class boundaryPointConstructor {
  private:
    // nada
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
    void set_point(point_t* p, boundaryType_t b);
    void delete_existing_point(point_t* p);
    void init_quader();
    void init_chopped_quader(point_t start, int devider);
    void init_quader(point_t start);
    void init_quader(point_t start,vector_t size);
    void init_sliding_lid();
    void init_chopped_sliding_lid(point_t start, int devider);
    void init_quader_side_chopped(point_t start, int chopsize);
    void init_sliding_lid_side_chopped(point_t start, int chopsize);
    void init_sliding_lid_inner(point_t start, point_t cntiues, vector_t inner_size);
    void delete_structures();
    // visualize
    void visualize_2D_boundary(int size);
    // total nodes
    long total_boundary_nodes();
};

#endif // NL_LATTICEBOLTZMANN_BOUNDARY_POINT_GENERATOR_H
