//
// Created by workstation on 16.11.22.
//

#ifndef MY_LATTICE_BOLTZMANN_H
#define MY_LATTICE_BOLTZMANN_H

// includes
#include "types.h"
#include "node.h"
// for completeness
extern matrix_t velocity_set;
// functions
typedef array_t (*equilibrium_function) (node* node);
typedef void (*streaming_function) (node* node); // could also do boundary conditions there?!
typedef void (*collision_function) (node* node);
typedef void (*calculate_macro_values) (node* node);



class simulation {
  private:
    std::vector<node*> nodes;
    int dimensions, channels;
    int size_x,size_y; // prob best to put them into arrays
    int limit_x, limit_y;
    //
    nodeIdentifier_t determine_node_type(int pox, int poy) const;
    void determine_neighbours();
    bool check_still_in_sim_space(array_t position) const;
  public:
    node* search_neighbour_node_body(node* hunter, array_t prey);
    node* search_neighbour_node_boundary(node* hunter);
    simulation() = default;
    void init(int six, int siy);
    void run();
    void get_data();
    std::vector<node*> access();
};

#endif // MY_LATTICE_BOLTZMANN_H
