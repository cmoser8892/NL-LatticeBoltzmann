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
    EXPECT_EQ(calculate_truncation_force(&c,&u,&f),0);
    c = {0,0};
    u = {1,-1};
    f = {1,1};
    EXPECT_EQ(calculate_truncation_force(&c,&u,&f),0);
    c = {1,0};
    u = {0,0};
    f = {1,1};
    EXPECT_EQ(calculate_truncation_force(&c,&u,&f),6);
    c = {1,0};
    u = {1,1};
    f = {1,1};
    EXPECT_EQ(calculate_truncation_force(&c,&u,&f),9);
    c = {1,1};
    u = {1,1};
    f = {1,1};
    EXPECT_EQ(calculate_truncation_force(&c,&u,&f),42);
}

TEST(FuntionalTest, conical_delta) {
    EXPECT_EQ(conical_delta(2,3),0);
    EXPECT_EQ(conical_delta(1,1),1);
}

// todo look up book boy Wolf Gladrow on forcing term in LB cap 5
// todo try out the guo term also described in viggen 6.14
// main equation is 5.2.9.
//  there are some strange cases still left -> investigate
// todo negative numbers in pkh?! do i even want that