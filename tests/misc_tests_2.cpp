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
    EXPECT_EQ(calculate_truncation_force(&c,&u,&f),3);
    c = {1,0};
    u = {1,1};
    f = {1,1};
    EXPECT_EQ(calculate_truncation_force(&c,&u,&f),6);
    c = {1,1};
    u = {1,1};
    f = {1,1};
    EXPECT_EQ(calculate_truncation_force(&c,&u,&f),36);
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

TEST(ForceTest, rotation_force_zero) {
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

TEST(ForceTest, zentrifugal_force) {
    point_t origin = {0,0};
    point_t canvas_size = {50,50};
    point_t test_point = {50,0};
    vector_t f;
    double omega = 1.00;
    goaForce test(origin,canvas_size,omega);
    test.calculate_F_rotation(0,0,&test_point);
    f = test.return_force_alpha();
    // std::cout << f << std::endl;
    // first component is just the centrifugal force
    // in the direction outboard dirction of the rotation
    // f_zentrifugal = w*w*r
    EXPECT_EQ(f.norm(),50);
    test_point = {50,50};
    test.calculate_F_rotation(0,0,&test_point);
    f = test.return_force_alpha();
    EXPECT_EQ(f.norm(),test_point.norm());
    omega = 2;
    test_point = {50,50};
    test.set_omega(omega);
    test.calculate_F_rotation(0,0,&test_point);
    f = test.return_force_alpha();
    EXPECT_EQ(f.norm(),4*test_point.norm());
}

TEST(ForceTest, coriolis_force) {
    point_t origin = {0,0};
    point_t canvas_size = {50,50};
    point_t test_point = {50,0};
    vector_t f;
    double omega = 1.00;
    goaForce test(origin,canvas_size,omega);
    test.calculate_F_rotation(1,0,&test_point);
    f = test.return_force_alpha();
    EXPECT_EQ(f.x(),0);
    EXPECT_EQ(f.y(),52);
    omega = 2.00;
    test.set_omega(omega);
    test_point =  {0,0};
    test.calculate_F_rotation(1,0,&test_point);
    f = test.return_force_alpha();
    EXPECT_EQ(f.x(),0);
    EXPECT_EQ(f.y(),4);
    omega = 2.00;
    test.set_omega(omega);
    test_point =  {0,0};
    test.calculate_F_rotation(0,1,&test_point);
    f = test.return_force_alpha();
    EXPECT_EQ(f.x(),-4);
    EXPECT_EQ(f.y(),0);
}

TEST(ForceTest, correct_truncation_terms) {
    // test against the c style implementation with the same inputs
    // the c function is assumed to be correct (tested previously)
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
        // std::cout << i << std::endl;
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
        // std::cout << i << std::endl;
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
        // std::cout << i << std::endl;
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

TEST(FunctionalTest, force_macro_calculation) {
    // simple value correction
    // setup arrays
    array_t forces;
    forces.resize(CHANNELS);
    forces.setZero();
    array_t values;
    values.resize(CHANNELS);
    values.setOnes();
    // sim class setup
    forcedSimulation fsim(nullptr, nullptr, nullptr);
    // we calculate some values
    {
        auto [rho,ux,uy] = fsim.test_calcualte_macro(&values,&forces);
        EXPECT_EQ(rho,9);
        EXPECT_EQ(ux, 0);
        EXPECT_EQ(uy, 0);
    }
    // add a force
    forces.setOnes();
    {
        auto [rho,ux,uy] = fsim.test_calcualte_macro(&values,&forces);
        EXPECT_EQ(rho,9+4.5);
        EXPECT_EQ(ux, 0);
        EXPECT_EQ(uy, 0);
    }
    {
        forces(0) = 0;
        forces(1) = 0;
        forces(2) = 0;
        forces(3) = 0;
        forces(4) = 0;
        auto [rho,ux,uy] = fsim.test_calcualte_macro(&values,&forces);
        EXPECT_EQ(rho,9+2);
        EXPECT_EQ(ux, 0);
        EXPECT_EQ(uy, 0);
    }
    {
        forces(2) = 1;
        forces(3) = 1;
        auto [rho,ux,uy] = fsim.test_calcualte_macro(&values,&forces);
        EXPECT_EQ(rho,9+3);
        EXPECT_EQ(ux, 0);
        EXPECT_EQ(uy, 0);
    }
    {
        forces.setOnes();
        forces(3) = 0;
        auto [rho,ux,uy] = fsim.test_calcualte_macro(&values,&forces);
        EXPECT_EQ(rho,9+4);
        EXPECT_EQ(ux, 0.5/rho);
        EXPECT_EQ(uy, 0.5/rho);
    }
}

TEST(FunctionalTest, force_term_first_second_moment) {
    // tests out weather or not the velocity moments of the truncation of the force term is calculated
    // correctly
    point_t origin = {0,0};
    point_t canvas_size = {50,50};
    double omega = 0.00;
    goaForce test(origin,canvas_size,omega);
    vector_t f = {12,8};
    vector_t v = {0,0};
    // test values
    double first_moment = 0;
    vector_t second_moment = {0,0};
    // set force
    test.set_force_alpha(f);
    test.set_velocity(v);
    test.calculate_F_i();
    // we calculate the 1st and 2nd moment of the force
    // 1st
    for(int i = 0; i < CHANNELS; ++i) {
        first_moment += test.force_channels[i]* weights(i);
    }
    // 2nd
    for(int i = 0; i < CHANNELS; ++i) {
        vector_t v_set = velocity_set.col(i);
        second_moment += v_set * (test.force_channels[i] * weights(i));
    }
    f = test.return_force_alpha();
    EXPECT_NEAR(second_moment.x(), f.x(),1e-10);
    EXPECT_NEAR(second_moment.y(), f.y(),1e-10);
    // std::cout << second_moment << std::endl;
    EXPECT_NEAR(first_moment, 0,1e-7);
}

TEST(FunctionalTest, equilibrium_moments) {
    double rho = 3;
    double ux = 5;
    double uy = 10;
    array_t f_eq = equilibrium_2d(ux, uy, rho);
    //
    double first_moment = 0;
    for(auto f: f_eq) {
        first_moment += f;
    }
    EXPECT_NEAR(first_moment,rho,1e-10);
    vector_t second_moment = {0,0};
    vector_t velocity_channel;
    for(int i = 0; i < CHANNELS; ++i) {
        velocity_channel = velocity_set.col(i);
        second_moment += velocity_channel * f_eq(i);
    }
    EXPECT_NEAR(second_moment(0)/rho,ux,1e-10);
    EXPECT_NEAR(second_moment(1)/rho,uy,1e-10);
}

// todo look up book boy Wolf Gladrow on forcing term in LB cap 5
// todo try out the guo term also described in viggen 6.14
// main equation is 5.2.9.
//  there are some strange cases still left -> investigate
// todo negative numbers in pkh?! do i even want that