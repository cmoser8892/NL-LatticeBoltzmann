
#ifndef NL_LATTICEBOLTZMANN_NODEGENERATOR_H
#define NL_LATTICEBOLTZMANN_NODEGENERATOR_H

/**
 * class to generate/hold all the node points and connect them
 * second ingredient to the init
 * takes the boundary points as an ingredient
 */
#include "types.h"
#include "boundary_point_generator.h"

typedef struct nodePoint {
    handle_t handle;
    array_t position;
    nodeIdentifier_t type;
    std::vector<toLinks_t*> links;
    boundaryType_t boundary; // ignored if node type is wet
}nodePoint_t;

/**
 * In the CG class we discussed memory arrangements for NL data,
 * a sort of Z shape seems most appropriate at least try to order them close in memory
 */
class nodeGenerator {
  private:
    boundaryPointConstructor* points;
    bool redo = false;
    void write_data_to_file();
    bool read_data_from_file();
    void read_back_switch_case();
    int determine_correct_channel();
    void determine_neighbors();
    void linear_generation();
    bool check_other_boundary_hit(boundaryPoint_t* p, point_t &check_point);
  public:
    std::vector<nodePoint_t*> node_infos;
    explicit nodeGenerator(boundaryPointConstructor* p);
    void init();
};

#endif // NL_LATTICEBOLTZMANN_NODEGENERATOR_H