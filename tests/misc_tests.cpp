// Include a library file to make sure proper includes are set
#include <gtest/gtest.h>
#include "lattice_boltzmann.h"

TEST(MiscTest, TestingWhatever) {
    int size_x = 5, size_y = 5;
    simulation test;
    test.init(size_x,size_y);
    std::cout << "init done" << std::endl;
    auto access = test.access();
}

TEST(InitTest, correctNodeAmount) {
    int size_x = 5, size_y = 5;
    simulation test;
    test.init(size_x,size_y);
    auto access = test.access();
    EXPECT_EQ(access.size(), size_x*size_y);
}

TEST(InitTest, boundaryNeighbors) {
    int size_x = 5, size_y = 5;
    simulation test;
    test.init(size_x,size_y);
    auto access = test.access();
    // test weather or not boundary neighours have the correct amount of nullptrs
    for(auto node : access) {
        if(node->node_type == BOUNDARY) {
            int counter = 0;
            for(auto n : node->neighbors) {
                if(n != nullptr) {
                    counter++;
                }
            }
            EXPECT_LE(counter,3);
        }
    }
}