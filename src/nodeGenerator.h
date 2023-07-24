#ifndef NL_LATTICEBOLTZMANN_NODEGENERATOR_H
#define NL_LATTICEBOLTZMANN_NODEGENERATOR_H

/**
 * class to generate/hold all the node points and connect them
 * second ingredient to the init
 * takes the boundary points as an ingredient
 */
#include "types.h"
#include "node.h"
#include "boundary_point_generator.h"
#include "straight.h"
#include "marker.h"

/**
 * Read back enum, to handle the lines.
 */
typedef enum readBack {
    HANDLE = 0, /**<  Handle first */
    TYPE, /**<  Dry Wet type */
    BOUNDARY, /**< Boundary type (Bounce Back, Moving, etc. ) */
    POSITION, /**<  Position of the node */
    LINKS, /**<   links to the node */
    ERROR = 255 /**< error is max uint8_t or char */
}readBack_t;

/**
 * In the CG class we discussed memory arrangements for NL data,
 * a sort of Z shape seems most appropriate at least try to order them close in memory
 */
class nodeGenerator {
  private:
    bool create_straight = true; /**< One variant creates the surfaces themself the other doesnt prevents double frees */
    boundaryPointConstructor* points = nullptr; /**<  Boundary point cloud pointer */
    std::string file_name = "stored_nodes_file"; /**< String to the name of the nodes */
    bool redo = true; /**<  Should we redo the node generator */
    bool save = false; /**<  Should we save */
    vector_t discovery_vector = {1, 0}; /**< Only relevant for linear creation */
    long size_canvas = 0;
    //
    void write_data_to_file(bool write);
    bool read_data_from_file();
    void read_back_switch_case(nodePoint_t* n, std::string& s, readBack_t *chop);
    void determine_neighbors();
    void linear_generation();
    bool check_other_boundary_hit(boundaryPoint_t* p, point_t &check_point);
    void check_nodes_inside();
    void check_nodes_ibm(double range);
    void remove_unwanted_nodes(handle_t* current);
    void add_boundary_nodes(handle_t* current);
    void reduce_boundary_neighborhood();
    void check_and_set_reduced_neighborhood(handle_t array_position, boundaryType_t b);
  public:
    std::vector<nodePoint_t*> node_infos; /**<  Node point link */
    std::vector<bool> to_be_removed; /**< Master control remove the node or not */
    straightGenerator* straight_surfaces = nullptr; /**< Straight surface pointer */
    markerIBM* markers = nullptr; /**< markers */
    explicit nodeGenerator(boundaryPointConstructor* p);
    explicit nodeGenerator(straightGenerator * s);
    ~nodeGenerator();
    void set_discovery_vector(vector_t set);
    void set_redo_save(bool r, bool s);
    void init();
    void init(unsigned int size);
    void init_fused(unsigned int size);
    void init_surface(unsigned int size, double range);
    void board_creation(unsigned int size);
    void delete_node_infos();
    void visualize_2D_nodes();
    void visualize_2D_nodes_labels(boundaryType_t t);
};

#endif // NL_LATTICEBOLTZMANN_NODEGENERATOR_H
