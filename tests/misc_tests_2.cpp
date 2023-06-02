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

TEST(ForceTest, circle_force) {
    double force = 0.0035;
    point_t dont_care = {0,0};
    point_t canvas_size = {50,50};
    goaForce test(dont_care,canvas_size,0);
    test.calculate_F_circle(&dont_care, force,0,0);
    // the norm of the force should be
    vector_t f = test.return_force_alpha();
    EXPECT_NEAR(f.norm(), force, 1e-5);
}

TEST(ForceTest, rotation_force) {
    point_t origin = {0,0};
    point_t canvas_size = {50,50};
    double omega = 0.00;
    goaForce test(origin,canvas_size,omega);
    test.calculate_F_rotation(0,0,&canvas_size);
    vector_t f = test.return_force_alpha();
    EXPECT_EQ(f.norm(),0);
    test.calculate_F_rotation(0,0,&origin);
    f = test.return_force_alpha();
    EXPECT_EQ(f.norm(),0);
}

TEST(ForceTest, functional_master_test_rotation_zero) {
    point_t origin = {0,0};
    point_t canvas_size = {50,50};
    double omega = 0.00;
    goaForce test(origin,canvas_size,omega);
    test.calculate_F_rotation(0,0,&canvas_size);
    test.calculate_F_i();
    for(auto channels : test.force_channels) {
        EXPECT_EQ(channels,0);
    }
}

TEST(ForceTest, correct_truncation_terms) {
    // test against the c style implementation with the same inputs
    point_t origin = {0,0};
    point_t canvas_size = {50,50};
    double omega = 0.00;
    goaForce test(origin,canvas_size,omega);
    vector_t f = {0,0};
    vector_t v = {0,0};
    // set force
    test.set_force_alpha(f);
    test.set_velocity(v);
    test.calculate_F_i();
    for(int i = 0; i < CHANNELS; ++i) {
        vector_t c = velocity_set.col(i);
        double check = calculate_truncation_force(&c,&v,&f);
        EXPECT_EQ(check,test.force_channels(i));
    }
    /// individual test
    f = {1,1};
    v = {1,1};
    test.set_force_alpha(f);
    test.set_velocity(v);
    test.calculate_F_i();
    for(int i = 0; i < CHANNELS; ++i) {
        std::cout << i << std::endl;
        vector_t c = velocity_set.col(i);
        double check = calculate_truncation_force(&c,&v,&f);
        EXPECT_NEAR(check,test.force_channels(i),1e-5);
    }
    /// individual test
    f = {2,3};
    v = {4,5};
    test.set_force_alpha(f);
    test.set_velocity(v);
    test.calculate_F_i();
    for(int i = 0; i < CHANNELS; ++i) {
        std::cout << i << std::endl;
        vector_t c = velocity_set.col(i);
        double check = calculate_truncation_force(&c,&v,&f);
        EXPECT_NEAR(check,test.force_channels(i),1e-5);
    }
    /// individual test
    f = {-1,-1};
    v = {-1,-1};
    test.set_force_alpha(f);
    test.set_velocity(v);
    test.calculate_F_i();
    for(int i = 0; i < CHANNELS; ++i) {
        std::cout << i << std::endl;
        vector_t c = velocity_set.col(i);
        double check = calculate_truncation_force(&c,&v,&f);
        EXPECT_NEAR(check,test.force_channels(i),1e-5);
    }
}

TEST(FuntionalTest, conical_delta) {
    EXPECT_EQ(conical_delta(2,3),0);
    EXPECT_EQ(conical_delta(1,1),1);
}

TEST(FunctionalTest, angles) {
    vector_t v1 = {1,0};
    vector_t v2 = {0,1};
    EXPECT_NEAR(calculate_angle(&v1, &v2), EIGEN_PI/2,1e-5);
    v2 = {1,1};
    EXPECT_NEAR(calculate_angle(&v1, &v2), EIGEN_PI/4,1e-5);
    v2 = {-1,0};
    EXPECT_NEAR(calculate_angle(&v1, &v2), EIGEN_PI,1e-5);
}

// todo look up book boy Wolf Gladrow on forcing term in LB cap 5
// todo try out the guo term also described in viggen 6.14
// main equation is 5.2.9.
//  there are some strange cases still left -> investigate
// todo negative numbers in pkh?! do i even want that