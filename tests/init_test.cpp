
#include <gtest/gtest.h>
#include "lattice_boltzmann.h"

TEST(InitTests, TestTemplate) {
    EXPECT_TRUE(true);
}

TEST(InitTests, TestNumberNeighbour) {
    // checks weather or not a body node has 8 neighbours
    int size_x = 5, size_y = 5;
    simulation test;
    test.init(size_x,size_y);
    auto access = test.access();
    for(auto node: access) {
        int counter = 0;
        for(auto i :node->neighbors) {
            if(i != nullptr)
                counter++;
        }
        std::cout << counter << std::endl;
        if(node->node_type == BODY)
            EXPECT_EQ(counter,8);
        // fails for boundary 6 neighbours?!
        if(node->node_type == BOUNDARY)
            EXPECT_EQ(counter,3);
    }
}