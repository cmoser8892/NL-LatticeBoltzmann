// Include a library file to make sure proper includes are set
#include <gtest/gtest.h>
#include "lattice_boltzmann.h"

TEST(MiscTest, TestingWhatever) {
    int size_x = 5, size_y = 5;
    simulation test;
    test.init(size_x,size_y);
    auto access = test.access();
    array_t s{0,0};
    auto whatever = test.search_neighbour_node(access.at(6),s);
}

TEST(InitTest, correctNodeAmount) {
    int size_x = 5, size_y = 5;
    simulation test;
    test.init(size_x,size_y);
    auto access = test.access();
    EXPECT_EQ(access.size(), size_x*size_y);
}
