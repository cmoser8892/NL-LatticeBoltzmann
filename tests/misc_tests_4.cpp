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
#include "drawn_image_surface.h"
#include <gtest/gtest.h>

/**
 * Tests some of the functionality of the ranging key lookup.
 * @test
 */
TEST(FunctionalTest, ranging_key_look) {
    //
    double range = 2;
    // Fill a ranging point key hash
    std::vector<point_t> points;
    handle_t runner = 0;
    rangingPointKeyHash rpkh;
    for(int x = 0; x < 5; ++x) {
        for(int y = 0; y < 5; ++y) {
            point_t temp = {x,y};
            points.push_back(temp);
            rpkh.fill_key(runner,temp);
            runner++;
        }
    }
    // do a ranging scan
    point_t middle  = {2,2};
    std::vector<handle_t> back = rpkh.ranging_key_translation(middle, range);
    EXPECT_EQ(back.size(),points.size());
    EXPECT_TRUE(rpkh.ranging_key_look_for_specific(middle,range,12));
    // add an additional point
    point_t additional = {3.2,3.2};
    rpkh.fill_key(runner,additional);
    EXPECT_TRUE(rpkh.ranging_key_look_for_specific(middle,2,runner));
}

/**
 * Tests the calculate  surface length function.
 * @test
 */
TEST(FunctionalTest, calcualte_surface_length) {
    straight_t input;
    straightGenerator sg;
    double side_length = 35; // with a distance of 0.75 we should get 80 markers
    point_t starter = {5,5};
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
    EXPECT_EQ(sg.calculate_total_surface_length(),4*35);
}

/**
 * Pretty much a bug case that i test against, going to the next surface seems bugged.
 * @test
 */
TEST(FunctionalTest, marker_go_next_overflow) {
    long canvas_size = 50;
    double marker_distance = 0.5;
    bool file_write = true;
    kernelType_t kernel = KERNEL_C;
    double ibm_distance = kernel_id_to_lattice_search(kernel);
    // Load the image
    straight_t input;
    straightGenerator sg;
    double side_length = 35; // with a distance of 0.75 we should get 80 markers
    point_t starter = {5,5};
    // we put in a quader
    input.point = starter;
    input.direction = {0,1};
    input.max_t = side_length;
    sg.add_surface(input);
    input.point += input.direction * side_length;
    vector_t v = {1,1};
    input.direction = v.normalized();
    input.max_t = 10;
    sg.add_surface(input);
    sg.surface_mass_center();
    nodeGenerator ng(&sg);
    ng.init_surface_return(canvas_size,ibm_distance,marker_distance);
    // surface calcs
    double surface_length = sg.calculate_total_surface_length();
    EXPECT_NEAR(surface_length,45,1e-5);
    // adjusted marker length so that an equal full marker length is filled
    double markers_fit_in= std::floor(surface_length/marker_distance);
    double add_on = std::fmod(surface_length,marker_distance);
    marker_distance += add_on/markers_fit_in;
    EXPECT_EQ(ng.markers->marker_points.size(),surface_length/marker_distance);
}