
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

typedef enum readBack {
    HANDLE = 0,
    TYPE,
    BOUNDARY,
    POSITION,
    LINKS,
    ERROR = 255
}readBack_t;
/**
 * In the CG class we discussed memory arrangements for NL data,
 * a sort of Z shape seems most appropriate at least try to order them close in memory
 */
class nodeGenerator {
  private:
    boundaryPointConstructor* points;
    std::string file_name = "stored_nodes_file";
    bool redo = true;
    bool save = false;
    vector_t discovery_vector = {1, 0};
    void write_data_to_file(bool write);
    bool read_data_from_file();
    void read_back_switch_case(nodePoint_t* n, std::string& s, readBack_t *chop);
    void determine_neighbors();
    void linear_generation();
    bool check_other_boundary_hit(boundaryPoint_t* p, point_t &check_point);
    void board_creation(unsigned int size);
  public:
    std::vector<nodePoint_t*> node_infos;
    explicit nodeGenerator(boundaryPointConstructor* p);
    void set_discovery_vector(vector_t set);
    void set_redo_save(bool r, bool s);
    void init();
    void init(unsigned int size);
};

#endif // NL_LATTICEBOLTZMANN_NODEGENERATOR_H
