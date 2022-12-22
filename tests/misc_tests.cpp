// Include a library file to make sure proper includes are set
#include "simulation.h"
#include "node.h"
#include "functions.h"
#include <gtest/gtest.h>

TEST(MiscTest, TestingWhatever) {
    EXPECT_TRUE(true);
}

TEST(FunctionalTest, correct_equilibrium) {
    // check weather or not the pen and computer agree
    handle_t h = 1;
    int dimension = 2;
    int channels = 9;
    array_t pos;
    pos.resize(3);
    pos << 1,2,4;
    auto n = new node(h,dimension,channels,pos);
    n->rho = 1;
    n->u.setZero();
    n->data = equilibrium(n);
    // should be standard init values
    EXPECT_NEAR(4.0/9,n->data(0),1e-10);
    EXPECT_NEAR(1.0/9,n->data(1),1e-10);
    EXPECT_NEAR(1.0/9,n->data(2),1e-10);
    EXPECT_NEAR(1.0/9,n->data(3),1e-10);
    EXPECT_NEAR(1.0/9,n->data(4),1e-10);
    EXPECT_NEAR(1.0/36,n->data(5),1e-10);
    EXPECT_NEAR(1.0/36,n->data(6),1e-10);
    EXPECT_NEAR(1.0/36,n->data(7),1e-10);
    EXPECT_NEAR(1.0/36,n->data(8),1e-10);
    // other values
    n->rho = 1;
    n->u.setOnes();
    n->data = equilibrium(n);
    EXPECT_NEAR(-8.0/9,n->data(0),1e-10);
    EXPECT_NEAR(11.0/18,n->data(1),1e-10);
    EXPECT_NEAR(11.0/18,n->data(2),1e-10);
    EXPECT_NEAR(-1.0/18,n->data(3),1e-10);
    EXPECT_NEAR(-1.0/18,n->data(4),1e-10);
    EXPECT_NEAR(22.0/36,n->data(5),1e-10);
    EXPECT_NEAR(-2.0/36,n->data(6),1e-10);
    EXPECT_NEAR(10.0/36,n->data(7),1e-10);
    EXPECT_NEAR(-2.0/36,n->data(8),1e-10);
    // -1 velocities
    n->rho = 1;
    n->u.setOnes();
    n->u *= -1;
    n->data = equilibrium(n);
    EXPECT_NEAR(-8.0/9,n->data(0),1e-10);
    EXPECT_NEAR(-1.0/18,n->data(1),1e-10);
    EXPECT_NEAR(-1.0/18,n->data(2),1e-10);
    EXPECT_NEAR(11.0/18,n->data(3),1e-10);
    EXPECT_NEAR(11.0/18,n->data(4),1e-10);
    EXPECT_NEAR(10.0/36,n->data(5),1e-10);
    EXPECT_NEAR(-2.0/36,n->data(6),1e-10);
    EXPECT_NEAR(22.0/36,n->data(7),1e-10);
    EXPECT_NEAR(-2.0/36,n->data(8),1e-10);
}

TEST(FunctionalTest,correct_macro) {
    // values should just read back
    handle_t h = 1;
    int dimension = 2;
    int channels = 9;
    array_t pos;
    pos.resize(3);
    pos << 1,2,4;
    auto n = new node(h,dimension,channels,pos);
    n->rho = 1;
    n->u.setZero();
    n->data = equilibrium(n);
    macro(n);
    EXPECT_NEAR(n->rho, 1,1e-10);
    EXPECT_NEAR(n->u(0), 0,1e-10);
    EXPECT_NEAR(n->u(1), 0,1e-10);
    // next
    n->rho = 1;
    n->u.setOnes();
    n->data = equilibrium(n);
    macro(n);
    EXPECT_NEAR(n->rho, 1, 1e-10);
    EXPECT_NEAR(n->u(0), 1, 1e-10);
    EXPECT_NEAR(n->u(1), 1, 1e-10);
    // next
    n->rho = 1;
    n->u.setOnes();
    n->u *= 10;
    n->data = equilibrium(n);
    macro(n);
    EXPECT_NEAR(n->rho, 1, 1e-10);
    EXPECT_NEAR(n->u(0), 10, 1e-10);
    EXPECT_NEAR(n->u(1), 10, 1e-10);
}

TEST(StreamTest, one_D_streaming_channel_one) {
    // for the streaming test the last two points have to be ignored (aka DRY nodes)
    int size = 20;
    point_t start = {0,0};
    point_t end = {size,0};
    point_t sim_area = {size,1};
    boundaryPointConstructor boundaries(sim_area);
    boundaries.set_point(&start,NO_BOUNDARY);
    boundaries.set_point(&end  ,NO_BOUNDARY);
    simulation sm(&boundaries);
    sm.init();
    // zero the data
    for(auto node : sm.nodes) {
        node->data.setZero();
    }
    //manual streaming channel 1 start at handle 1 at 1,0 to handle 4 at 4,0
    sm.nodes.at(0)->data(1) = 1;
    for(int i = 1; i < size; ++i) {
        // go check
        for(auto node : sm.nodes) {
            int expect = 0;
            if(node->handle == i) {
                expect = 1;
            }
            EXPECT_EQ(expect, node->data(1));
        }
        // execute two steps
        sm.streaming_step_1();
        sm.streaming_step_2();
    }
}

TEST(StreamTest, one_D_streaming_channel_three) {
    int size = 20;
    point_t start = {0,0};
    point_t end = {size,0};
    point_t sim_area = {size,1};
    boundaryPointConstructor boundaries(sim_area);
    boundaries.set_point(&start,NO_BOUNDARY);
    boundaries.set_point(&end  ,NO_BOUNDARY);
    simulation sm(&boundaries);
    sm.init();
    // zero the data
    for(auto node : sm.nodes) {
        node->data.setZero();
    }
    //manual streaming channel 1
    // the -3 is wierd, but if u think about it last tow nodes
    // have handle 20 and 19 in array 19 and 18 so last on is 17 in array
    sm.nodes.at(size-3)->data(3) = 1;
    std::cout << sm.nodes.at(size-3)->handle << std::endl;
    // aka 20 nodes in total numbered 1 to 20 last two a boundary so we ignore those
    for(int i = size-2; i > 0; --i) {
        // go check
        for(auto node : sm.nodes) {
            // std::cout << node->neighbors.size() << std::endl;
            int expect = 0;
            if(node->handle == i) {
                expect = 1;
            }
            EXPECT_EQ(expect, node->data(3));
        }
        // execute two steps
        sm.streaming_step_1();
        sm.streaming_step_2();
    }
}

TEST(BounceBackTesting, South) {
    // we need to generate two nodes and only really consider the bb
}
