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

/**
 * The init for ibm simulations is tested here, we want the right amount of nodes and so on.
 * @test
 */
TEST(IbmTest, init_ibm) {
    // set up a straight generator to do the full init in the new ibm class as the current one is just a mess
    straightGenerator sg;
    point_t starter = {3.5,3.5};
    straight_t input;
    long canvas_size = 13;
    double side_length = 5; // with a distance of 0.75 we should get 80 markers
    double ibm_distance = 2;
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
    sg.surface_mass_center();
    // generator
    nodeGenerator ng(&sg);
    ng.init_surface(canvas_size,ibm_distance);
    ng.visualize_2D_nodes();
    EXPECT_EQ(ng.node_infos.size(),81);
    // force
    point_t dk = {0,0};
    vector_t sizes = {canvas_size,canvas_size};
    goaForce rot(dk,sizes,1e-5);
    // sim init
    ibmSimulation sim(&ng,&rot, ng.markers, sizes);
    sim.init();
    EXPECT_EQ(sim.nodes.size(), 81);
    EXPECT_EQ(sim.markers.size(), ng.markers->marker_points.size());
    // check ibm nodes
    long ibm_numbers = 0;
    for(auto n : sim.nodes) {
        if(n->boundary_type == IBM) {
            ibm_numbers++;
        }
    }
    // EXPECT_EQ(ibm_numbers, 4*(8+6+4+2)); // only the center is not set
    handle_t mh = 0;
    for(auto m : sim.markers) {
        EXPECT_EQ(++mh,m->handle);
    }
}

/**
 * Tests if all the nodes have the correct amount of links
 * @test
 * @attention sim nodes got to have 8 links, links that have less in the node info get self pointers to the node
 */
TEST(IbmTest, nodes_placement_links) {
    // init a structure
    point_t starter = {3,3};
    straight_t input;
    straightGenerator sg;
    long canvas_size = 13;
    double side_length = 6; // with a distance of 0.75 we should get 80 markers
    double ibm_distance = 2;
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
    sg.surface_mass_center();
    nodeGenerator ng(&sg);
    ng.init_surface(canvas_size,ibm_distance);
    ng.visualize_2D_nodes();
    simulation_parameters params;
    params.ibm_range = ibm_distance;
    point_t dk = {0,0};
    vector_t sizes = {canvas_size,canvas_size};
    goaForce rot(dk,sizes,1e-3);
    ibmSimulation sim(&ng, &rot,ng.markers,sizes);
    sim.init();
    // check right init
    int ibm_nodes = 0;
    int pure_nodes = 0;
    int undefined = 0;
    for(auto n : sim.nodes) {
        if(n->boundary_type == IBM) {
            ibm_nodes++;
        }
        else if(n->boundary_type == NO_BOUNDARY) {
            pure_nodes++;
        }
        else {
            undefined++;
        }
    }
    // note for equal on point it is times 5 otherwise not really lol
    EXPECT_EQ(ibm_nodes,4*(side_length*5));
    EXPECT_EQ(pure_nodes, sim.nodes.size() - ibm_nodes);
    EXPECT_EQ(undefined, 0);
    // check links in the sim everyone should be a 8
    array_t link_sizes;
    link_sizes.setZero(9);
    for(auto n : sim.nodes) {
        link_sizes(n->neighbors.size())++;
    }
    EXPECT_EQ(link_sizes(0),0);
    EXPECT_EQ(link_sizes(1),0);
    EXPECT_EQ(link_sizes(2),0);
    EXPECT_EQ(link_sizes(3),0);
    EXPECT_EQ(link_sizes(4),0);
    EXPECT_EQ(link_sizes(5),0);
    EXPECT_EQ(link_sizes(6),0);
    EXPECT_EQ(link_sizes(7),0);
    EXPECT_EQ(link_sizes(8),sim.nodes.size());
    // checks sizes in the node generator
    link_sizes.setZero();
    for(auto n : ng.node_infos) {
        link_sizes(n->links.size())++;
    }
    EXPECT_EQ(link_sizes(0),0);
    EXPECT_EQ(link_sizes(1),0);
    EXPECT_EQ(link_sizes(2),0);
    EXPECT_EQ(link_sizes(3),4);
    EXPECT_EQ(link_sizes(4),0);
    EXPECT_EQ(link_sizes(5),4*9);
    EXPECT_EQ(link_sizes(6),0);
    EXPECT_EQ(link_sizes(7),0);
    EXPECT_EQ(link_sizes(8),ng.node_infos.size()-40);
}

/**
 * Tests the marker movement while doing simulation.
 * @test
 * @note bug somewhere seems to emit from the lower left corner
 */
TEST(IbmTest, marker_movement_around) {
    // init a structure
    point_t starter = {3,3};
    straight_t input;
    straightGenerator sg;
    long canvas_size = 13;
    double side_length = 6; // with a distance of 0.75 we should get 80 markers
    double ibm_distance = 2;
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
    sg.surface_mass_center();
    nodeGenerator ng(&sg);
    ng.init_surface(canvas_size,ibm_distance);
    ng.visualize_2D_nodes();
    simulation_parameters params;
    params.ibm_range = ibm_distance; // has to be corrected
    params.k = 1;
    point_t dk = {0,0};
    vector_t sizes = {canvas_size,canvas_size};
    goaForce rot(dk,sizes,1e-3);
    ibmSimulation sim(&ng, &rot,ng.markers,sizes);
    sim.init();
    // check right init
    int ibm_nodes = 0;
    int pure_nodes = 0;
    int undefined = 0;
    for(auto n : sim.nodes) {
        if(n->boundary_type == IBM) {
            ibm_nodes++;
        }
        else if(n->boundary_type == NO_BOUNDARY) {
            pure_nodes++;
        }
        else {
            undefined++;
        }
    }
    // note for equal on point it is times 5 otherwise not really lol
    EXPECT_EQ(ibm_nodes,4*(side_length*5));
    EXPECT_EQ(pure_nodes, sim.nodes.size() - ibm_nodes);
    EXPECT_EQ(undefined, 0);
    // check for the right marker movement
    int steps = 0;
    sim.run(0);
    // check where the markers are
    for(auto m : sim.markers) {
        std::cout << m->original_position.x() << " ," << m->original_position.y() << std::endl;
        EXPECT_NEAR(m->original_position.x(), m->position.x(), 1e-1);
        EXPECT_NEAR(m->original_position.y(), m->position.y(), 1e-1);
    }
}

/**
 * Tests the individual movement of a marker.
 * @test
 */
TEST(IbmTest, marker_movement_individual_set) {
    // init a structure
    point_t starter = {3,3};
    straight_t input;
    straightGenerator sg;
    long canvas_size = 13;
    double side_length = 6; // with a distance of 0.75 we should get 80 markers
    double ibm_distance = 2;
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
    sg.surface_mass_center();
    nodeGenerator ng(&sg);
    ng.init_surface(canvas_size,ibm_distance);
    // ng.visualize_2D_nodes();
    simulation_parameters params;
    params.ibm_range = ibm_distance;
    params.k = 10;
    point_t dk = {0,0};
    vector_t sizes = {canvas_size,canvas_size};
    goaForce rot(dk,sizes,1e-3);
    ibmSimulation sim(&ng, &rot,ng.markers,sizes);
    sim.init();
    // note we abuse the init here just to get markers and i want to test the integration of
    vector_t velocity = {1,1};
    for(auto m : sim.markers) {
        m->velocity = velocity;
    }
    sim.test_propagate_markers();
    for(auto m : sim.markers) {
        EXPECT_EQ(velocity+m->original_position,m->position);
    }
    for(auto m : sim.markers) {
        m->velocity = velocity;
    }
    sim.test_propagate_markers();
    for(auto m : sim.markers) {
        EXPECT_EQ(2*velocity+m->original_position,m->position);
    }
}

/**
 * Tests correct kernel function and calcualtes the integral under the function.
 * @test
 */
TEST(IbmTest, kernels) {
    double kernel_stencil = 1;
    point_t origin = {1,1};
    point_t next = {0,0};
    vector_t r = next - origin;
    //
    double factorised_kernel = (kernel_3(r.x(),kernel_stencil)*kernel_3(r.y(),kernel_stencil))/(kernel_stencil*kernel_stencil);
    double un_factorised_kernel = kernel_3(r.norm(),kernel_stencil);
    EXPECT_NE(un_factorised_kernel,factorised_kernel);
    double func = d_kernel_32(&r,kernel_stencil);
    EXPECT_EQ(func, factorised_kernel);
    // integration test
    // set up x and y
    int total = 240;
    double begin = -2.5;
    double end = 2.5;
    double range = end - begin;
    double step = range/total;
    array_t x;
    x.resize(total);
    array_t y;
    y.resize(total);
    for(int i = 0; i < total; ++i) {
        x[i] = begin;
        y[i] = begin;
        begin += step;
    }
    // calculate z
    double integral = 0;
    flowfield_t z;
    z.resize(total,total);
    for(int i = 0; i < total; ++i) {
        for(int j = 0; j < total; ++j) {
            vector_t dr = {x(i),y(j)};
            z(i,j) =  d_kernel_32(&dr,kernel_stencil);
            integral += z(i,j) * step * step;
        }
    }
    EXPECT_NEAR(integral, 1,1e-5);
    // check that the other one is not working
    integral = 0;
    z.setZero();
    for(int i = 0; i < total; ++i) {
        for(int j = 0; j < total; ++j) {
            vector_t dr = {x(i),y(j)};
            z(i,j) = kernel_3(dr.norm(),kernel_stencil);
            integral += z(i,j) * step * step;
        }
    }
    EXPECT_NE(1,integral); // about 1.87 so complete nonsense
}

TEST(IbmTest, self_stream) {
    // tests the self stream in nodes that have less than 8 links
}

TEST(IbmTest, force_strength) {
    // tests for correct values in the force implementation
}


// todo there are some strange cases still left -> investigate
// todo refactor boundary point generators look for bumps and so on -> ranging pgk is a great tool her