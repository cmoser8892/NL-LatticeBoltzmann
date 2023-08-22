
/*
 * vallgrind call
valgrind --tool=callgrind --dump-instr=yes (p)
 */

#include "helper_functions.h"
#include "drawn_image_surface.h"
#include "nodeGenerator.h"

#include <iostream>
#include <fstream>

/**
 * 00 main, testing stuff.
 * @note good rule of thumb is all 4 contours need to be selected
 * @return
 */
int main() {
    //
    std::ofstream out;
    out.open("out.txt");
    // board creation for some z maps
    straightGenerator* sg = nullptr;
    nodeGenerator ng(sg);
    ng.board_creation(48);
    int i = 0;
    for(auto n : ng.node_infos) {
        ++i;
        out << "[" << n->position.x() << ", " <<n->position.y() << "], ";
        if((i % 10) == 0) {
            out << std::endl;
        }
    }
    out.close();
}
