// Include a library file to make sure proper includes are set
#include "simulation.h"
#include "node.h"
#include "functions.h"
#include "helper_functions.h"
#include "neighborhood.h"
#include "image_converter.h"
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


// todo investiagate the odd 0 pass in the intersection tests write out a full test for that
// todo there are some strange cases still left -> in vestigate
// todo negative numbers in pkh?! do i even want that