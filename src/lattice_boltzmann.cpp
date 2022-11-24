#include "lattice_boltzmann.h"
#include "functions.h"
#include "helper_functions.h"
#include <iostream>

/// nodes
// constructor should be the only thing needed
node::node(int dimensions, int channels, array_t pos, node_identifier_t type,int ar_pos) {
    node_type = type;
    rho = 1;
    data.resize(channels);
    copy.resize(channels);
    u.resize(dimensions);
    position = pos;
    array_position = ar_pos;
}

/// simulation run class
node_identifier_t simulation::determine_node_type(int pox, int poy) {
    node_identifier_t return_value = BODY;
    if( pox == 0 || poy == 0 || pox == limit_x || poy == limit_y) {
        return_value = BOUNDARY;
    }
    return return_value;
}

void simulation::determine_neighbours() {
    // neighbours 1 to channels
    for (auto node : nodes) {
        // i gives the channel number
        for(int i = 1; i < node->data.size();++i) {
            // this is prob the most lazy implementation ever
            array_t search;
            search.resize(dimensions);
            if( node->node_type == BODY) {
                search = node->position + velocity_set(i);
            }
            if (node->node_type == BOUNDARY) {
                search = node->position - velocity_set(i);
            }
            node->neighbors.push_back(search_neighbour_node(node,search));
        }
    }
}

node* simulation::search_neighbour_node(node *hunter, array_t prey) {
    // basic old search function fallthrough?!
    node* return_node = nullptr;
    if (check_still_in_sim_space(prey))
    {
        // get own postion
        array_t search_position = hunter->position;
        // calculate a counter form the hunter to the prey
        array_t bullet = prey - search_position;
        int counter = int(round(bullet.x() + bullet.y()*size_x));
        // set up an iterator for quick reference
        auto reference_iter = nodes.begin()+hunter->array_position;
        // add the counter calculated beforehand should be the right point
        reference_iter += counter; // crashes if outside of the normal domain
        auto m = reference_iter.operator*();
        if(compare_arrays(m->position, prey)) {
            // only include body nodes
            if(m->node_type == BODY)
                return_node = m;
        }
        else {
            std::cout << "Couldnt find neighbor node, try again" << std::endl;
            // brute force search function as a fallback
            for(auto s : nodes) {
                if((s->position(0) == prey(0)) && (s->position(1) == prey(1))){
                    return_node = s;;
                    break;
                }
            }
        }
    }
    return return_node;
}

bool simulation::check_still_in_sim_space(array_t position) {
    bool return_value = true;
    if(position.x() < 0)
        return_value = false;
    if(position.y() < 0)
        return_value = false;
    if(position.x() > limit_x)
        return_value = false;
    if(position.y() > limit_y)
        return_value = false;
    return return_value;
}

void simulation::init(int six, int siy) {
    size_x = six; size_y = siy;
    limit_x = six-1; limit_y = siy -1;
    dimensions = 2;
    channels = 9;
    // create nodes
    int counter = 0;
    for( int y = 0; y < siy; ++y) {
        for(int x = 0; x < six; ++x) {
            array_t p = {{double(x)},{double(y)}};
            node_identifier_t type = determine_node_type(x,y);
            node* n = new node(dimensions,channels, p, type,counter);
            nodes.push_back(n);
            counter++;
        }
    }
    // determine neighbours
    determine_neighbours();
}

void simulation::run() {
    for(auto node: nodes) {
        // streaming includes bounce back in the most basic form (no u wall thou todo)
        streaming_step1(node);
        streaming_step2(node);
        macro(node);
        collision(node);
    }
}

std::vector<node*> simulation::access() {
    return nodes;
}


