#include "node_generator.h"
#include <iostream>

node_generator::node_generator(boundaryPointConstructor *p) {
    points = p;
}

void node_generator::linear_generation() {
    int handle_counter = 1;
    // go throu the boundary points starting at a b point and go throu while still discovering new ones

    // lastlly add the boundary points
    for(auto p : points->boundary_points) {
        std::cout << p << std::endl;
        auto n = new nodePoint_t;
        n->handle = handle_counter;
        n->position = p->point;
        n->type = DRY;
        n->boundary = p->type;
    }
}
