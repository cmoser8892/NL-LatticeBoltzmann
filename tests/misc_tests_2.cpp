// Include a library file to make sure proper includes are set
#include "simulation.h"
#include "node.h"
#include "functions.h"
#include "helper_functions.h"
#include "neighborhood.h"
#include "image_converter.h"
#include <gtest/gtest.h>

// the misc file got to laggy apparently i should stop at 2000 lines
TEST(ForceTest, rotation_eigen) {
    Eigen::Rotation2D<double> rot;
    rot.angle() = EIGEN_PI/2;
    vector_t test;
    test << 1, 0;
    test = rot * test;
    EXPECT_NEAR(test.y(),1,1e-10);
    EXPECT_NEAR(test.x(),0,1e-10);
}

TEST(ForceTest, truncation_force) {
    // setup the vectors
    vector_t c = {0,0};
    vector_t u = {0,0};
    vector_t f = {0,0};
    EXPECT_EQ(calculate_truncation_force(c,u,f),0);
    c = {0,0};
    u = {1,-1};
    f = {1,1};
    EXPECT_EQ(calculate_truncation_force(c,u,f),0);
    c = {1,0};
    u = {0,0};
    f = {1,1};
    EXPECT_EQ(calculate_truncation_force(c,u,f),6);
    c = {1,0};
    u = {1,1};
    f = {1,1};
    EXPECT_EQ(calculate_truncation_force(c,u,f),9);
    c = {1,1};
    u = {1,1};
    f = {1,1};
    EXPECT_EQ(calculate_truncation_force(c,u,f),42);
}

TEST(ForceTest, rotating_zeros) {
    point_t origin = {0,0};
    point_t middle = {5,5};
    point_t point = middle;
    rotatingForce rf(origin,middle,0,0);
    rf.precalculate(0,0,&point);
    EXPECT_EQ(rf.return_force_alpha().x(),0);
    EXPECT_EQ(rf.return_force_alpha().y(),0);
}

TEST(ForceTest, vector_force) {
    // a bit of a simple and non-memory efficient way to implement different forces
    int total_size = 160;
    int swicht = 10;
    double magnitude = 0.1;
    std::vector<vector_t> f = circular_force_generation(total_size,swicht,magnitude);
    EXPECT_EQ(f.size(),total_size);
    for(auto v : f) {
        EXPECT_EQ(v.norm(),magnitude);
    }
}

TEST(ForceTest, circular_force) {
    // test for the circular class implementation
    int swicht = 10;
    double magnitude = 0.1;
    circularForce force(swicht,magnitude);
    for(int i = 0; i < 9; ++i) {
        EXPECT_EQ(magnitude,force.return_current_x());
        EXPECT_EQ(magnitude,force.return_next_x());
        EXPECT_EQ(0, force.return_current_y());
        EXPECT_EQ(0, force.return_next_y());
        force.increment();
    }
    EXPECT_EQ(magnitude,force.return_current_x());
    EXPECT_EQ(0,force.return_next_x());
    EXPECT_EQ(0, force.return_current_y());
    EXPECT_EQ(magnitude, force.return_next_y());
    force.increment();
    EXPECT_EQ(0 ,force.return_current_x());
    EXPECT_EQ(0 ,force.return_next_x());
    EXPECT_EQ(magnitude, force.return_current_y());
    EXPECT_EQ(magnitude, force.return_next_y());
}

// todo look up book boy Wolf Gladrow on forcing term in LB cap 5
// todo try out the guo term also described in viggen 6.14
// main equation is 5.2.9.
// todo update fused init to also work with boundary limits,
//  there are some strange cases still left -> investigate
// todo read about std::atomic may be useful
// todo negative numbers in pkh