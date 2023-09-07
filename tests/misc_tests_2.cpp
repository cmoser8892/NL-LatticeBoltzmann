// Include a library file to make sure proper includes are set
#include "simulation.h"
#include "node.h"
#include "functions.h"
#include "helper_functions.h"
#include "neighborhood.h"
#include "image_converter.h"
#include "marker.h"
#include "forces.h"
#include "lbm_simulation.h"
#include <gtest/gtest.h>
// the misc file got to laggy apparently i should stop at 2000 lines

/**
 * Checks out how ration with eigen work.
 * @test
 */
TEST(ForceTest, rotation_eigen) {
    Eigen::Rotation2D<double> rot;
    rot.angle() = EIGEN_PI/2;
    vector_t test;
    test << 1, 0;
    test = rot * test;
    EXPECT_NEAR(test.y(),1,1e-10);
    EXPECT_NEAR(test.x(),0,1e-10);
}

/**
 * Tests weather or not the implementation of the truncation force is correct.
 * @test
 * @see calculate_truncation_force()
 */
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
/**
 * Tests if the force calculation is correct.
 * @test
 * @see goaForce::calculate_F_circle()
 */
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

/**
 * Tests if the rotation force is 0 if all the ingredients are 0.
 * @test
 * @see goaForce::calculate_F_rotation()
 */
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

/**
 * Tests that the trunication terms are 0 for a force of 0.
 * @test
 * @see goaForce::calculate_F_rotation()
 * @see goaForce::calculate_F_i()
 */
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

/**
 * Tests against expected values for the zentrifugal force.
 * @test
 * @see goaForce::calculate_F_rotation()
 */
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

/**
 * Tests weather or not the coriolis force is calculated against expected values.
 * @test
 * @see goaForce::calculate_F_rotation()
 */
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

/**
 * Tests the correct truncation function terms agaist expected values.
 * @test
 * @see calculate_truncation_force()
 */
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
    // individual test
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
    // individual test
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
    // individual test
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

/**
 * Tests the conical delta function.
 * @test
 * @see conical_delta()
 */
TEST(FuntionalTest, conical_delta) {
    EXPECT_EQ(conical_delta(2,3),0);
    EXPECT_EQ(conical_delta(1,1),1);
}

/**
 * Tests teh calculate_angle function against expected values.
 * @test
 * @see calculate_angle()
 */
TEST(FunctionalTest, angles) {
    vector_t v1 = {1,0};
    vector_t v2 = {0,1};
    EXPECT_NEAR(calculate_angle(&v1, &v2), EIGEN_PI/2,1e-5);
    v2 = {1,1};
    EXPECT_NEAR(calculate_angle(&v1, &v2), EIGEN_PI/4,1e-5);
    v2 = {-1,0};
    EXPECT_NEAR(calculate_angle(&v1, &v2), EIGEN_PI,1e-5);
}

/**
 * Test the force sim method for calculation of macro values against expected values.
 * @test
 * @see forcedSimulation::test_calcualte_macro()
 */
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
        EXPECT_EQ(rho,9+0.5);
        EXPECT_EQ(ux, 0);
        EXPECT_EQ(uy, 0);
    }
}

/**
 * Checks weather or not the first and second moment of the truncation term of the force are correct.
 * @test
 * @see goaForce::calculate_F_i()
 */
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

/**
 * Tests weather or not the moments of the equilibrium function are correct.
 * @test
 * @see equilibrium_2d()
 */
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

/**
 * Tests the modified macro calculation in a force simulation against another way of implementaiton..
 * @test
 * @see forcedSimulation::test_calcualte_macro()
 * @see calculate_macro_population()
 */
TEST(FunctionalTest, macro_tests) {
    // stand in for the populations
    array_t values;
    values.resize(CHANNELS);
    values.setOnes();
    // sim class setup
    point_t origin = {0,0};
    point_t canvas_size = {50,50};
    double omega = 0.00;
    goaForce test(origin,canvas_size,omega);
    forcedSimulation fsim(nullptr, nullptr, &test);
    // force and velocity
    vector_t f = {0,0};
    vector_t v = {0,0};
    // we calculate some values
    {
        // we have no effect from the force
        test.set_force_alpha(f);
        test.set_velocity(v);
        test.calculate_F_i();
        auto [rho,ux,uy] = fsim.test_calcualte_macro(&values,&test.force_channels);
        auto [test_rho,test_ux, test_uy] = calculate_macro_population(&values);
        EXPECT_EQ(rho,test_rho);
        EXPECT_EQ(ux, test_ux);
        EXPECT_EQ(uy, test_uy);
    }
    {
        // we have no effect from the force
        f = {3,0};
        v = {5,6};
        test.set_force_alpha(f);
        test.set_velocity(v);
        test.calculate_F_i();
        auto [rho,ux,uy] = fsim.test_calcualte_macro(&values,&test.force_channels);
        auto [test_rho,test_ux, test_uy] = calculate_macro_population(&values);
        // test values
        double first_moment = 0;
        vector_t second_moment = {0,0};
        // 1st
        for(int i = 0; i < CHANNELS; ++i) {
            first_moment += test.force_channels[i]* weights(i);
        }
        // 2nd
        for(int i = 0; i < CHANNELS; ++i) {
            vector_t v_set = velocity_set.col(i);
            second_moment += v_set * (test.force_channels[i] * weights(i));
        }
        EXPECT_EQ(rho,test_rho);
        // in debug i can use eq in rel i have to uses near cause optimization
        EXPECT_NEAR(ux, test_ux + 0.5*f.x()/rho,1e-15);
        EXPECT_NEAR(uy, test_uy + 0.5*f.y()/rho,1e-15);
    }
}

/**
 * Tests the function to find multiple keys in case there are more.
 * @note there are two variants with a coordinate and a point (its just the same function with some boilerplate)
 * @test
 * @see pointKeyHash::multi_key_translation()
 */
TEST(FunctionalTest, multiple_pkh_entries) {
    pointKeyHash pkh;
    int number = 10;
    double var = 0.1; // how much the points are appart
    // create a vector of points
    std::vector<point_t> storage_points;
    point_t p = {0,0};
    vector_t v = {1,1};
    for(int i = 0; i < number; ++i) {
        p += i*var*v;
        pkh.fill_key(i,p);
        storage_points.push_back(p);
    }
    // test function in its intended use
    for(int i = 0; i < number; ++i) {
        std::vector<handle_t> handles = pkh.multi_key_translation(storage_points[i]);
        int counter = 0;
        for(auto h : handles) {
            if(storage_points[i] == storage_points[h]) {
                counter++;
            }
        }
        EXPECT_EQ(counter,1);
    }
}

/**
 * Tests the ranging functions basic functionality
 * @test
 * @see pointKeyHash::ranging_key_translation()
 */
TEST(FunctionalTest, handle_ranges_basic) {
    pointKeyHash pkh;
    int number = 30;
    double var = 0.1; // how much the points are appart
    point_t middle = {1.5,1.5};
    // create a vector of points
    std::vector<point_t> storage_points;
    point_t p = {0,0};
    vector_t v = {1,1};
    for(int i = 0; i < number; ++i) {
        pkh.fill_key(i,p);
        storage_points.push_back(p);
        p += var*v;
    }
    /// test w
    std::vector<handle_t> test;
    test = pkh.ranging_key_translation(middle,1);
    EXPECT_EQ(test.size(),number);
}

/**
 * Tests with more cells filled.
 * @test
 * @see pointKeyHash::ranging_key_translation()
 */
TEST(FunctionalTest, every_cell_filled) {
    pointKeyHash pkh;
    int number = 30;
    double var = 0.1; // how much the points are appart
    point_t middle = {1.5,1.5};
    // create a vector of points
    std::vector<point_t> storage_points;
    point_t p = {0,0};
    vector_t v = {1,1};
    for(int i = 0; i < number; ++i) {
        // std::cout << p << std::endl;
        pkh.fill_key(i,p);
        storage_points.push_back(p);
        p += var*v;
    }
    p = {0,3};
    v = {1,-1};
    for(int i = 30; i < number*2; ++i) {
        // std::cout << p << std::endl;
        pkh.fill_key(i,p);
        storage_points.push_back(p);
        p += var*v;
    }
    std::vector<handle_t> test;
    test = pkh.ranging_key_translation(middle,1);
    EXPECT_EQ(test.size(),number*2-1); // first one is not part
}

/**
 * Tests the correct composition of cells (weather or not all get actually touched).
 * @test
 * @see pointKeyHash::ranging_key_translation()
 */
TEST(FunctionalTest, one_point_in_cell) {
    pointKeyHash pkh;
    int number = 9;
    point_t middle = {12.5,12.5};
    // go over vector
    std::vector<point_t> storage_points;
    for(int i = 0; i < number; ++i) {
        vector_t vel_set = velocity_set.col(i);
        point_t c = middle + vel_set;
        pkh.fill_key(i,c);
        storage_points.push_back(c);
    }
    // find the guys
    std::vector<handle_t> test;
    test = pkh.ranging_key_translation(middle,1);
    EXPECT_EQ(number,test.size());
}

/**
 * Tests the ranging pkh functionality over a bigger range.
 * @test
 * @see pointKeyHash::ranging_key_translation()
 */
TEST(FunctionalTest, bigger_pkh_ranging) {
    pointKeyHash pkh;
    int number = 49;
    point_t middle = {3.5,3.5};
    // go over vector
    std::vector<point_t> storage_points;
    int i = 1;
    for(int x = 0; x < 7; ++x) {
        for(int y = 0; y < 7; ++y) {
            point_t c = {x,y};
            pkh.fill_key(i,c);
            storage_points.push_back(c);
            ++i;
        }
    }
    // find the guys
    std::vector<handle_t> test;
    test = pkh.ranging_key_translation(middle,3);
    EXPECT_EQ(number,test.size());
}

/**
 * Tests the right functionality of a ranging scan over a boarder.
 * @note Pkh expects only positive values
 * @see pointKeyHash::ranging_key_translation()
 */
TEST(FunctionalTest, over_pkh_zero_limit) {
    pointKeyHash pkh;
    int number = 4;
    point_t middle = {0.5,0.5};
    // go over vector
    std::vector<point_t> storage_points;
    int i = 1;
    for(int x = 0; x < 2; ++x) {
        for(int y = 0; y < 2; ++y) {
            point_t c = {x,y};
            pkh.fill_key(i,c);
            storage_points.push_back(c);
            ++i;
        }
    }
    // find the guys
    std::vector<handle_t> test;
    test = pkh.ranging_key_translation(middle,1);
    EXPECT_EQ(number,test.size());
}

/**
 * The standard langragian use case where we have a number of points and an lagragin point searching.
 * @test
 * @see pointKeyHash::ranging_key_translation
 */
TEST(FunctionalTest, lagragian_use_case) {
    pointKeyHash pkh;
    int number = 40;
    // go over vector
    std::vector<point_t> storage_points;
    int i = 1;
    for(int x = 0; x < 40; ++x) {
        for(int y = 0; y < 40; ++y) {
            point_t c = {x,y};
            pkh.fill_key(i,c);
            storage_points.push_back(c);
            ++i;
        }
    }
    // we search in a standard 2 wide search (for one of the kernel functions)
    point_t lagragian_point = {12.2,13.2};
    int range = 2;
    std::vector<handle_t> test = pkh.ranging_key_translation(lagragian_point, range);
    EXPECT_EQ(test.size(),25);
}

/**
 * Tests what happens with a ranging search for a zero search length.
 * @test
 * @see pointKeyHash::ranging_key_translation
 */
TEST(FunctionalTest, zero_ranging_pgk) {
    pointKeyHash pkh;
    int number = 4;
    point_t middle = {0.5,0.5};
    // go over vector
    std::vector<point_t> storage_points;
    int i = 1;
    for(int x = 0; x < number; ++x) {
        for(int y = 0; y < number; ++y) {
            point_t c = {x,y};
            pkh.fill_key(i,c);
            storage_points.push_back(c);
            ++i;
        }
    }
    std::vector<handle_t> test = pkh.ranging_key_translation(middle, 0);
    // we should find one (we search the origin cell)
    EXPECT_EQ(test.size(),1);
}

/**
 * Tests the envisioned use case for a marker force disipoation.
 * @test
 * @see rangingPointKeyHash::ranging_key_translation()
 */
TEST(FunctionalTest, simple_lagra_use_case) {
    rangingPointKeyHash rpkh;
    // go over vector
    std::vector<point_t> storage_points;
    int i = 1;
    for(int x = 0; x < 5; ++x) {
        for(int y = 0; y < 5; ++y) {
            point_t c = {x,y};
            rpkh.fill_key(i,c);
            storage_points.push_back(c);
            ++i;
        }
    }
    // we search in a standard 2 wide search (for one of the kernel functions)
    point_t lagragian_point = {2.5,2.5};
    int range = 2;
    std::vector<handle_t> test = rpkh.ranging_key_translation(lagragian_point, range);
    EXPECT_EQ(16,test.size());
}

/**
 * Just a check weather or not adding coordinates works.
 * @test
 * @see add_coordinates()
 */
TEST(FunctionalTest, coord_add) {
    coordinate_t a;
    coordinate_t b;
    a.x = 2;
    a.y = 3;
    b.x = 4;
    b.y = 5;
    coordinate_t c = add_coordinates(a,b);
    EXPECT_EQ(c.x, 6);
    EXPECT_EQ(c.y, 8);
}

/**
 * Basic tests for the kernel 1 function.
 * @test
 * @see kernel_A()
 */
TEST(FunctionalTest, kernel_A) {
    double range = 0;
    double range_x = 1;
    range = 0;
    EXPECT_EQ(kernel_A(range),1);
    range = 1;
    EXPECT_EQ(kernel_A(range),0);
    range = -1;
    EXPECT_EQ(kernel_A(range),0);
    range = 3;
    EXPECT_EQ(kernel_A(range),0);
    range = -4;
    EXPECT_EQ(kernel_A(range),0);
}

/**
 * Basic test for the kernel 2 function.
 * @test
 * @see kernel_B()
 */
TEST(FunctionalTest, kernel_B) {
    double range = 0;
    double range_x = 1;
    range = 0;
    EXPECT_EQ(kernel_B(range),2.0/3);
    range = -1.5;
    EXPECT_EQ(kernel_B(range),0);
    range = 1.5;
    EXPECT_EQ(kernel_B(range),0);
    range = -3;
    EXPECT_EQ(kernel_B(range),0);
    range = 3;
    EXPECT_EQ(kernel_B(range),0);
}

/**
 * Basic test for the kernel_C function.
 * @test
 * @see kernel_C()
 */
TEST(FunctionalTest, kernel_C) {
    double range = 0;
    double range_x = 1;
    range = 0;
    EXPECT_EQ(kernel_C(range),4.0/8);
    range = 2;
    EXPECT_EQ(kernel_C(range),0);
    range = -2;
    EXPECT_EQ(kernel_C(range),0);
    range = 5;
    EXPECT_EQ(kernel_C(range),0);
    range = -5;
    EXPECT_EQ(kernel_C(range),0);
}

/**
 * Follow up test where there is no actual fluid node and everything is quite close.
 * @test
 * @note no fluid notes are actually found bit of a blind test
 */
TEST(FunctionalTest, no_fluid_boundary) {
    unsigned int size = 2;
    unsigned int sub_size = 2;
    point_t c = {size,size}; // size or the canvas
    point_t p = {sub_size,sub_size}; // size of the  quader
    boundaryPointConstructor boundaries(c);
    // boundaries.init_sliding_lid_side_chopped({20,10},30);
    boundaries.init_quader({0,0},p);
    // construct straights
    straightGenerator s(&boundaries);
    s.init();
    nodeGenerator n(&boundaries);
    n.init_fused(size);
    EXPECT_EQ(n.node_infos.size(), 0);
}

/**
 * We have on straight on which we distribute markers.
 * @test
 * @see markerIBM::distribute_markers()
 */
TEST(FunctionalTest, distribute_markers_base) {
    straight_t input;
    straightGenerator sg;
    // we put in a quader
    input.point = {1,1};
    input.direction = {0,1};
    input.max_t = 6;
    input.type = IBM;
    sg.add_surface(input);
    markerIBM m(&sg);
    m.distribute_markers();
    EXPECT_EQ(m.marker_points.size(),8);
}

/**
 * We make some quaders and test weather or not the code does what it should.
 * @test
 * @see markerIBM::distribute_markers()
 */
TEST(FunctionalTest, distribute_markers_quader) {
    double marker_dist = 0.75;
    point_t starter = {1.4,1.4};
    for(int i = 96; i < 120; ++i) {
        straightGenerator sg;
        // we want sth without any residual and % only works on int
        double fit_in = (double) i / marker_dist;
        double floored_fit_in = std::floor(fit_in);
        if(std::fmod(fit_in,floored_fit_in) == 0) {
            // side length
            double side_length = (double) i / 4;
            double expected_value = fit_in;
            straight_t input;
            // we put in a quader
            input.point = starter;
            input.direction = {0,1};
            input.max_t = side_length;
            input.type = IBM;
            sg.add_surface(input);
            input.point += input.direction * side_length;
            input.direction = {1,0};
            input.max_t = side_length;
            input.type = IBM;
            sg.add_surface(input);
            input.point += input.direction * side_length;
            input.direction = {0,-1};
            input.max_t = side_length;
            input.type = IBM;
            sg.add_surface(input);
            input.point += input.direction * side_length;
            input.direction = {-1,0};
            input.max_t = side_length;
            input.type = IBM;
            sg.add_surface(input);
            // init the markers
            markerIBM mibm(&sg);
            mibm.distribute_markers();
            EXPECT_EQ(mibm.marker_points.size(),expected_value);
        }
    }
}
/**
 * We construct a quader that can not fit the 0.75 standard distance.
 * @test
 */
TEST(FunctionalTest, odd_quader) {
    point_t starter = {1,1};
    straight_t input;
    straightGenerator sg;
    double side_length = 5;
    // we put in a quader
    input.point = starter;
    input.direction = {0,1};
    input.max_t = side_length;
    input.type = IBM;
    sg.add_surface(input);
    input.point += input.direction * side_length;
    input.direction = {1,0};
    input.max_t = side_length;
    input.type = IBM;
    sg.add_surface(input);
    input.point += input.direction * side_length;
    input.direction = {0,-1};
    input.max_t = side_length;
    input.type = IBM;
    sg.add_surface(input);
    input.point += input.direction * side_length;
    input.direction = {-1,0};
    input.max_t = side_length;
    input.type = IBM;
    sg.add_surface(input);
    // tests how many will get generatored
    markerIBM mibm(&sg);
    mibm.distribute_markers();
    EXPECT_EQ(mibm.marker_points.size(),26);
    EXPECT_GT(mibm.return_marker_distance(),0.75);
}

/**
 * Tests/Explores the functionality when we actually have two structures.
 * @test
 * @attention based on the result of this test one can assume that surfaces can be structure blind.
 * @note still have to be in order otherwise results wont work (it might still be a good idea to partition thou -> parrallel)
 */
TEST(FunctionalTest, two_structure) {
    point_t starter = {1,1};
    straight_t input;
    straightGenerator sg;
    double side_length = 15; // with a distance of 0.75 we should get 80 markers
    // we put in a quader
    input.point = starter;
    input.direction = {0,1};
    input.max_t = side_length;
    input.type = IBM;
    sg.add_surface(input);
    input.point += input.direction * side_length;
    input.direction = {1,0};
    input.max_t = side_length;
    input.type = IBM;
    sg.add_surface(input);
    input.point += input.direction * side_length;
    input.direction = {0,-1};
    input.max_t = side_length;
    input.type = IBM;
    sg.add_surface(input);
    input.point += input.direction * side_length;
    input.direction = {-1,0};
    input.max_t = side_length;
    input.type = IBM;
    sg.add_surface(input);
    // another one
    starter = {4,4};
    side_length = 5;
    // we put in a quader
    input.point = starter;
    input.direction = {0,1};
    input.max_t = side_length;
    sg.add_surface(input);
    input.point += input.direction * side_length;
    input.direction = {1,0};
    input.max_t = side_length;
    sg.add_surface(input);
    input.point += input.direction * side_length;
    input.direction = {0,-1};
    input.max_t = side_length;
    sg.add_surface(input);
    input.point += input.direction * side_length;
    input.direction = {-1,0};
    input.max_t = side_length;
    sg.add_surface(input);
    // main thing i do not want to see is a whole when we change the surface
    // easiest way for this to work is to set up a ranging pkh and look for the markers
    // number of different cases still but should be doable
    markerIBM mibm(&sg);
    mibm.distribute_markers();
    EXPECT_EQ(mibm.marker_points.size(),80+26); // we got the right amount
    array_t on_surface;
    on_surface.resize(sg.surfaces.size());
    on_surface.setZero();
    double dk;
    std::vector<int> surface_structure_1;
    std::vector<int> surface_structure_2;
    for(int i = 0; i < mibm.marker_points.size(); ++i) {
        auto m = mibm.marker_points[i];
        for(int i = 0; i < sg.surfaces.size(); ++i) {
            auto s = sg.surfaces[i];
            if(point_on_straight(s,m,&dk)) {
                on_surface[i]++;
                if(i < 4) {
                    surface_structure_1.push_back(i);
                }
                else {
                    surface_structure_2.push_back(i);
                }
            }
        }
    }
    // check out the spacing on each individual surface
    for(auto n : surface_structure_1) {
        int next = (n+1) % (int)surface_structure_1.size();
        double distance = (*mibm.marker_points[n] - *mibm.marker_points[next]).norm();
        EXPECT_NEAR(mibm.return_marker_distance(), distance,1e-8);
    }
    for(auto n : surface_structure_2) {
        int next = (n+1) % (int)surface_structure_2.size();
        double distance = (*mibm.marker_points[n] - *mibm.marker_points[next]).norm();
        EXPECT_NEAR(mibm.return_marker_distance(), distance,1e-8);
    }
}

/**
 * Point on straight test.
 * @test
 * @attention negatives are not tested cause they dont appear in the code
 * @see point_on_straight()
 */
TEST(FunctionalTest, points_on_straights) {
    straight_t test;
    test.point = {0,0};
    test.direction = {1,1};
    test.max_t = 5;
    // positives
    double dk = 0;
    point_t test_point = {0,0};
    EXPECT_TRUE(point_on_straight(&test,&test_point, &dk));
    test_point = {3,3};
    EXPECT_TRUE(point_on_straight(&test,&test_point, &dk));
    test_point = {5,5};
    EXPECT_TRUE(point_on_straight(&test,&test_point, &dk));
    // negatives
    test_point = {6,6};
    EXPECT_TRUE(!point_on_straight(&test,&test_point, &dk));
    EXPECT_EQ(dk, 1);
    // different straight
    test.direction = {2,1};
    test.direction.normalize();
    test_point = {2,1};
    EXPECT_TRUE(point_on_straight(&test,&test_point,&dk));
    test_point = {4,2};
    EXPECT_TRUE(point_on_straight(&test,&test_point,&dk));
    test_point = {6,3};
    EXPECT_TRUE(!point_on_straight(&test,&test_point,&dk));
    test.direction = {0,1};
    test_point = {0,2};
    EXPECT_TRUE(point_on_straight(&test,&test_point,&dk));
    test_point = {1,2};
    EXPECT_TRUE(!point_on_straight(&test,&test_point,&dk));
    test.direction = {1,0};
    test_point = {2,0};
    EXPECT_TRUE(point_on_straight(&test,&test_point,&dk));
    // understand norm in eigen
    vector_t n = test.direction.normalized();
    EXPECT_NEAR(n.norm(),1,1e-8);
}

/**
 * Not really a tests just a idk.
 * @test
 */
TEST(FunctionalTest, SpringForces) {
    double k = 0.01;
    double d = 0.75;
    double dx = 1;
    point_t original_point = {0,0};
    point_t current_point = {-1,-1};
    vector_t norm = current_point - original_point;
    vector_t force = -k * d/dx * norm;
    EXPECT_TRUE(true);
}

/**
 * tests the truncation force.
 * @test
 */
TEST(FunctionalTest, trunication_force) {
    vector_t f = {3,5};
    vector_t u = {2,9};
    point_t dk = {0,0};
    goaForce g(dk,dk,0);
    g.set_force_alpha(f);
    g.set_velocity(u);
    g.calculate_F_i();
    array_t test = calculate_truncation_array(&f,&u);
    for(int i = 0; i < test.size(); ++i) {
        test[i] = g.force_channels[i];
    }
    EXPECT_TRUE(true);
}
