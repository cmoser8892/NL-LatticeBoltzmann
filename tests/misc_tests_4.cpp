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
#include "helper_functions.h"
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
    // Load the image
    straight_t input;
    straightGenerator sg;
    double side_length = 35; // with a distance of 0.75 we should get 80 markers
    point_t starter = {5,5};
    // we put in a quader
    input.point = starter;
    input.direction = {0,1};
    input.max_t = side_length;
    input.type = IBM;
    sg.add_surface(input);
    input.point += input.direction * side_length;
    vector_t v = {1,1};
    input.direction = v.normalized();
    input.max_t = 10;
    input.type = IBM;
    sg.add_surface(input);
    sg.surface_mass_center();
    nodeGenerator ng(&sg);
    ng.init_surface_return(canvas_size,kernel,marker_distance);
    // surface calcs
    double surface_length = sg.calculate_total_surface_length();
    EXPECT_NEAR(surface_length,45,1e-5);
    // adjusted marker length so that an equal full marker length is filled
    double markers_fit_in= std::floor(surface_length/marker_distance);
    double add_on = std::fmod(surface_length,marker_distance);
    marker_distance += add_on/markers_fit_in;
    EXPECT_EQ(ng.markers->marker_points.size(),surface_length/marker_distance);
}

/**
 * Tests out how we can identify the edges.
 * @test
 */
TEST(FunctionalTest, open_boundaries_id_test) {
    // node generator variables
    long canvas_size = 50;
    vector_t draw_size = {50-1,50-1};
    double marker_distance = 0.5;
    bool file_write = true;
    kernelType_t kernel = KERNEL_C;
    // Load the image
    auto test_image = get_base_path();
    test_image.append("tests");
    test_image.append("test_images");
    test_image.append("black_bars.png");
    // call the drawer
    surfaceDrawer s(test_image);
    std::vector<int> sel = {0,1};
    s.run_non_connecting(sel, false);
    s.close_open_surface(draw_size);
    s.surface_storage.surface_mass_center();
    nodeGenerator ng(&s.surface_storage);
    ng.init_surface_return(canvas_size,kernel,marker_distance);
    // check if we get the right amount of surfaces
    EXPECT_EQ(s.surface_storage.surfaces.size(),4);
}

/**
 * Compares the oldest and the newest collision implementation
 * @todo there might be an mistake 1/relaxation somewhere, it is really confusing relaxation tau and omega are used indistinguishable
 * @test
 */
TEST(FunctionalTest, collision_againt_test) {
    point_t dk = {1,1};
    node n(1,2,9,dk,NO_BOUNDARY);
    ibmSimulation test(nullptr,nullptr,nullptr,dk);
    double rho = 1;
    double ux = 1;
    double uy = 1;
    double relaxation = 0.5;
    test.parameters.relaxation = 1/relaxation;
    // add a node
    fNode fn(1,9,NO_BOUNDARY);
    fn.position = dk;
    // set the node variables
    n.rho = rho;
    n.u(0) = ux;
    n.u(1) = uy;
    // function calls
    test.test_collision(&fn,rho,ux,uy);
    collision(&n,relaxation);
    // compare values
    for(int i = 0; i < CHANNELS; ++i) {
        EXPECT_NEAR(fn.populations(i),n.population_even(i),1e-10);
    }
}

/**
 * Preliminary node test when working with periodics.
 * @test
 */
TEST(FunctionalTest, tube_no_markers) {
    // node generator variables
    long canvas_size = 150;
    vector_t draw_size = {canvas_size-1,canvas_size-1};
    double marker_distance = 0.5;
    bool file_write = true;
    kernelType_t kernel = KERNEL_C;
    // introduce a surface
    straight_t input;
    straightGenerator lines;
    // set the input and introduce
    input.point = {0,27.1};
    input.direction = {1,0};
    input.max_t = 149;
    input.type = IBM;
    lines.add_surface(input);
    input.point = {149,82.8};
    input.direction = {-1,0};
    input.max_t = 149;
    input.type = IBM;
    lines.add_surface(input);
    // introduce the containing lines will have to be manually found !
    input.point = {0,27.1};
    input.direction = {0,1};
    input.max_t = 82.8-27.1;
    input.type = PERIODIC;
    lines.add_surface(input);
    input.point = {149,27.1};
    input.direction = {0,1};
    input.max_t = 82.8-27.1;
    input.type = PERIODIC;
    lines.add_surface(input);
    // setup the node generator
    nodeGenerator ng(&lines);
    ng.init_surface_return(canvas_size,KERNEL_C,marker_distance);
    // ng.visualize_2D_nodes();
    EXPECT_EQ(ng.node_infos.size(),150*63);
    //
    long counter_periodic = 0;
    long counter_ibm_p = 0;
    long counter_none = 0;
    long counter_default = 0;
    for(auto n : ng.node_infos) {
        if(n->type == PERIODIC_CONNECT) {
            ++counter_periodic;
            EXPECT_EQ(n->links.size(),8);
            switch(n->boundary) {
            case INIT_NONE: {
                counter_none++;
                break;
            }
            case IBM_OUTER: {
                counter_ibm_p++;
                break;
            }
            default:
                counter_default++;
            }
        }
    }
    // checks
    EXPECT_EQ(2*63,counter_periodic);
    EXPECT_EQ(4*8,counter_ibm_p);
    EXPECT_EQ(2*63-4*8,counter_none);
    EXPECT_EQ(counter_default,0);
}

/**
 * Method to place periodic markers on a surface.
 * @test
 */
TEST(FunctionalTest, periodidics_marker_placement) {
    straight_t input;
    input.point = {0,27.1};
    input.direction = {0,1};
    input.max_t = 82.8-27.1;
    input.type = PERIODIC;

    markerPoints periodic_markers(nullptr,1);
    periodic_markers.distribute_markers_periodic(&input,NO_BOUNDARY,KERNEL_A);
    EXPECT_EQ(periodic_markers.marker_points.size(),1 + 54);
    // check if empty and so on
    if(!periodic_markers.marker_points.empty()) {
        for (size_t i = 0; i < (int(periodic_markers.marker_points.size()) - 1); ++i) {
            auto current = periodic_markers.marker_points[i];
            auto next = periodic_markers.marker_points[i +1];
            // test the distance between them
            vector_t t = *next - *current;
            // checks the distance between the marker points should be 1
            EXPECT_NEAR(t.norm(),1,1e-10);
        }
    }
    // check the first point and last points
    point_t test = {0,28};
    EXPECT_EQ(test.y(),periodic_markers.marker_points[0]->y());
    test = {0,82};
    EXPECT_EQ(test.y(),periodic_markers.marker_points[periodic_markers.marker_points.size()-1]->y());
}

/**
 * Method to place periodic markers on a surface.
 * @test
 */
TEST(FunctionalTest, periodidics_marker_placement_ibm_usecase) {
    straight_t input;
    input.point = {0,27.1};
    input.direction = {0,1};
    input.max_t = 82.8-27.1;
    input.type = PERIODIC;
    // this and that
    markerPoints periodic_markers(nullptr,1);
    periodic_markers.distribute_markers_periodic(&input,IBM,KERNEL_C);
    EXPECT_EQ(periodic_markers.marker_points.size(),1 + 54 + 8);
    // check if empty and so on
    if(!periodic_markers.marker_points.empty()) {
        for (size_t i = 0; i < (int(periodic_markers.marker_points.size()) - 1); ++i) {
            auto current = periodic_markers.marker_points[i];
            auto next = periodic_markers.marker_points[i +1];
            // test the distance between them
            vector_t t = *next - *current;
            // checks the distance between the marker points should be 1
            EXPECT_NEAR(t.norm(),1,1e-10);
        }
    }
    // check the first point and last points
    point_t test = {0,24};
    EXPECT_EQ(test.y(),periodic_markers.marker_points[0]->y());
    test = {0,86};
    EXPECT_EQ(test.y(),periodic_markers.marker_points[periodic_markers.marker_points.size()-1]->y());
}

/**
 * Tests the check_plus_minus_90 function.
 * @tests
 */
TEST(FunctionalTest, ninty_angle_off) {
    vector_t ref = {-1,0};
    vector_t test = {0,0};
    // tests
    test = {0,0};
    EXPECT_TRUE(check_plus_minus_90(&test,&ref));
    test = {-1,0};
    EXPECT_TRUE(check_plus_minus_90(&test,&ref));
    test = {0,1};
    EXPECT_TRUE(check_plus_minus_90(&test,&ref));
    test = {0,-1};
    EXPECT_TRUE(check_plus_minus_90(&test,&ref));
    test = {1,0};
    EXPECT_TRUE(!check_plus_minus_90(&test,&ref));
    test = {1,1};
    EXPECT_TRUE(!check_plus_minus_90(&test,&ref));
    test = {1,-1};
    EXPECT_TRUE(!check_plus_minus_90(&test,&ref));
    test = {-1,-1};
    EXPECT_TRUE(check_plus_minus_90(&test,&ref));
    test = {-1,1};
    EXPECT_TRUE(check_plus_minus_90(&test,&ref));
    test = {-500000,0};
    EXPECT_TRUE(check_plus_minus_90(&test,&ref));
}

/**
 * Tests the recenter of a reference vector to a cardinal vector.
 * @test
 */
TEST(FunctionalTest, min_cardinal) {
    vector_t ref = {0,1.001};
    vector_t expected = {0,1};
    vector_t test = vector_to_cardinal(ref);
    // vectors and points are fundamentally the same object
    EXPECT_TRUE(compare_two_points(&expected,&test));
}

/**
 * Tests the vector to index in the velocity set function.
 * @test
 */
TEST(FunctionalTest, set_to_index) {
    vector_t test = {0,0};
    vector_t compare = {0,0};
    // tests good
    test = {0,0};
    EXPECT_EQ(index_of_velocity_set(test),0);
    test = {1,0};
    EXPECT_EQ(index_of_velocity_set(test),1);
    test = {0,1};
    EXPECT_EQ(index_of_velocity_set(test),2);
    test = {-1,0};
    EXPECT_EQ(index_of_velocity_set(test),3);
    test = {0,-1};
    EXPECT_EQ(index_of_velocity_set(test),4);
    test = {1,1};
    EXPECT_EQ(index_of_velocity_set(test),5);
    test = {-1,1};
    EXPECT_EQ(index_of_velocity_set(test),6);
    test = {-1,-1};
    EXPECT_EQ(index_of_velocity_set(test),7);
    test = {1,-1};
    EXPECT_EQ(index_of_velocity_set(test),8);
    // bad
    test = {3,0};
    EXPECT_EQ(index_of_velocity_set(test),-1);
}

/**
 * Tests an i swap expression.
 * @test
 */
TEST(FunctionalTest, little_things) {
    int i = 0;
    int k = (i+1)%2;
    EXPECT_EQ(k,1);
    i = 1;
    k = (i+1)%2;
    EXPECT_EQ(k,0);
}

TEST(FunctionalTest, inlet_to_outlet_test) {
    // node generator variables
    long canvas_size = 150;
    vector_t draw_size = {canvas_size-1,canvas_size-1};
    double marker_distance = 0.5;
    bool file_write = true;
    kernelType_t kernel = KERNEL_C;
    // introduce a surface
    straight_t input;
    straightGenerator lines;
    // set the input and introduce
    input.point = {0,27.1};
    input.direction = {1,0};
    input.max_t = 149;
    input.type = IBM;
    lines.add_surface(input);
    input.point = {149,82.8};
    input.direction = {-1,0};
    input.max_t = 149;
    input.type = IBM;
    lines.add_surface(input);
    // introduce the containing lines will have to be manually found !
    straight_t inlet;
    inlet.point = {0,27.1};
    inlet.direction = {0,1};
    inlet.max_t = 82.8-27.1;
    inlet.type = PERIODIC;
    lines.add_surface(inlet);
    straight_t outlet;
    outlet.point = {149,27.1};
    outlet.direction = {0,1};
    outlet.max_t = 82.8-27.1;
    outlet.type = PERIODIC;
    lines.add_surface(outlet);
    // setup the node generator
    nodeGenerator ng(&lines);
    ng.init_surface_return(canvas_size,KERNEL_C,marker_distance);
    // ng.visualize_2D_nodes();
    // check creation
    EXPECT_EQ(ng.node_infos.size(),150*63);
    // create the marker points
    markerPoints inlet_markers(nullptr,1);
    inlet_markers.distribute_markers_periodic(&inlet,IBM,KERNEL_C);
    markerPoints outlet_markers(nullptr,1);
    outlet_markers.distribute_markers_periodic(&outlet,IBM,KERNEL_C);
    // setup a sim
    vector_t sizes = {canvas_size,canvas_size};
    ibmSimulation tester(&ng, nullptr,ng.markers,sizes);
    tester.init();
    // get to the inlet and outlet
    // create a ranging pkh
    rangingPointKeyHash rpkh;
    handle_t current = 1;
    for(auto node : tester.nodes) {
        rpkh.fill_key(current,node->position);
        ++current;
    }
    std::vector<handle_t> inlet_handles = {};
    for(auto point : inlet_markers.marker_points) {
        std::vector<handle_t> temp = rpkh.ranging_key_translation(*point,0.9);
        // put into inlet handles vector
        if(temp.size() == 1) {
            inlet_handles.push_back(temp[0]);
        }
        else {
            // trigger error
            EXPECT_TRUE(false);
        }
    }
    // italicising pasta again
    std::vector<handle_t> outlet_handles = {};
    for(auto point : outlet_markers.marker_points) {
        std::vector<handle_t> temp = rpkh.ranging_key_translation(*point,0.9);
        // put into inlet handles vector
        if(temp.size() == 1) {
            outlet_handles.push_back(temp[0]);
        }
        else {
            // trigger error
            EXPECT_TRUE(false);
        }
    }
    // zero stuff
    for(auto n : tester.nodes) {
        n->populations.setZero();
    }
    // set values
    for(auto h: inlet_handles) {
        auto current_node = tester.nodes[h-1];
        // set population
        current_node->populations(3) = 3;
        current_node->populations(CHANNELS + 3) = 3;
    }
    for(auto h : outlet_handles) {
        auto current_node = tester.nodes[h-1];
        // set population
        current_node->populations(1) = 1;
        current_node->populations(CHANNELS + 1) = 1;
    }
    // do some streaming
    tester.test_streaming(0);
    int population_position = tester.offset_sim;
    // do the actual testing
    for(auto h: inlet_handles) {
        auto current_node = tester.nodes[h-1];
        // check
        EXPECT_EQ(current_node->populations(population_position + 1),1);
    }
    for(auto h: outlet_handles) {
        auto current_node = tester.nodes[h-1];
        // check
        EXPECT_EQ(current_node->populations(population_position + 3),3);
    }

}