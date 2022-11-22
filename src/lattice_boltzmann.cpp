#include "lattice_boltzmann.h"
#include "functions.h"
matrix_t velocity_set = {{0,1,0,-1,0 ,1,-1,-1, 1},
                         {0,0,1,0 ,-1,1, 1,-1,-1}};


/// nodes
// constructor should be the only thing needed
node::node(int dimensions, int channels, array_t pos, node_identifier_t type) {
    node_type = type;
    rho = 1;
    data.resize(channels);
    copy.resize(channels);
    u.resize(dimensions);
    position = pos;
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
            bool found = false;
            search.resize(dimensions);
            if( node->node_type == BODY) {
                search = node->position + velocity_set(i);
            }
            if (node->node_type == BOUNDARY) {
                search = node->position - velocity_set(i);
                // pox = 0   need 367 to 185
                // poy = 0   need 478 to 256
                // pox = lim need 158 to 376
                // poy = lim need 256 to 478
            }
            // search function kinda lazy i know
            for(auto s : nodes) {
                if((s->position(0) == search(0)) && (s->position(1) == search(1))){
                    node->neighbors.push_back(s);
                    found = true;
                }
            }
            if(found == false) {
                node->neighbors.push_back(nullptr);
            }
        }
    }
}

void simulation::init(int six, int siy) {
    size_x = six; size_y = siy;
    limit_x = six-1; limit_y = siy -1;
    dimensions = 2;
    channels = 9;
    // create nodes
    for( int x = 0; x < six; ++x) {
        for(int y = 0; y < siy; ++y) {
            array_t p = {{double(x)},{double(y)}};
            node_identifier_t type = determine_node_type(x,y);
            node* n = new node(dimensions,channels, p, type);
            nodes.push_back(n);
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


