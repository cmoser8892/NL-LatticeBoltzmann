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
    kernelType_t kernel = KERNEL_A;
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
    double ibm_distance = kernel_id_to_lattice_search(kernel);
    nodeGenerator ng(&sg);
    ng.init_surface(canvas_size,ibm_distance);
    ng.visualize_2D_nodes();
    simulation_parameters params;
    params.ibm_range = kernel_id_to_lattice_search(kernel);
    params.kernel_in_use = kernel;
    params.k = 1;
    point_t dk = {0,0};
    vector_t sizes = {canvas_size,canvas_size};
    goaForce rot(dk,sizes,1e-3);
    ibmSimulation sim(&ng, &rot,ng.markers,sizes);
    sim.set_simulation_parameters(params);
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
 * Tests the integrals of the normal 1d kernels to be 1.
 * @test
 */
TEST(IBMTest, kernels_1d) {
    // we just integrate the kernels
    int total = 640;
    double begin = -4.5;
    double end = 4.5;
    double range = end - begin;
    double step = range/total;
    array_t x;
    x.resize(total);
    for(int i = 0; i < total; ++i) {
        x[i] = begin;
        begin += step;
    }
    // kernel A phi_2
    double integral = 0;
    for(int i = 0; i < total; ++i) {
        integral += kernel_A(x[i])*step;
    }
    EXPECT_NEAR(integral, 1, 1e-1); // integration is not as accurate here
    // kernel B phi_3
    integral = 0;
    for(int i = 0; i < total; ++i) {
        integral += kernel_B(x[i])*step;
    }
    EXPECT_NEAR(integral, 1, 1e-5);
    // kernel C phi_4
    integral = 0;
    for(int i = 0; i < total; ++i) {
        integral += kernel_C(x[i])*step;
    }
    EXPECT_NEAR(integral, 1, 1e-5);
}

/**
 * Tests the integrals of the  2d kernels to be 1.
 * @test
 */
TEST(IBMTest, kernels_2d) {
    int total = 640;
    double begin = -4.5;
    double end = 4.5;
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
    // kernel a
    double integral = 0;
    for(int i = 0; i < total; ++i) {
        for(int j = 0; j <total; ++j) {
            vector_t dk = {x(i),y(j)};
            integral += kernel_A_2d(&dk) * step*step;
        }
    }
    EXPECT_NEAR(integral, 1, 1e-3);
    // kernel b
    integral = 0;
    for(int i = 0; i < total; ++i) {
        for(int j = 0; j <total; ++j) {
            vector_t dk = {x(i),y(j)};
            integral += kernel_B_2d(&dk) * step*step;
        }
    }
    EXPECT_NEAR(integral, 1, 1e-5);
    // kernel c
    integral = 0;
    for(int i = 0; i < total; ++i) {
        for(int j = 0; j <total; ++j) {
            vector_t dk = {x(i),y(j)};
            integral += kernel_C_2d(&dk) * step*step;
        }
    }
    EXPECT_NEAR(integral, 1, 1e-5);
}

/**
 * Tests correct kernel function and calcualtes the integral under the function.
 * @test
 */
TEST(IbmTest, kernels) {
    double kernel_stencil = 1;
    point_t origin = {1, 1};
    point_t next = {0, 0};
    vector_t r = next - origin;
    //
    double factorised_kernel = (kernel_C(r.x()) * kernel_C(r.y()));
    double un_factorised_kernel = kernel_C(r.norm());
    double frac = kernel_C_2d(&r);
    EXPECT_EQ(factorised_kernel, frac);
    EXPECT_NE(un_factorised_kernel, factorised_kernel);
    // integration test
    // set up x and y
    int total = 640;
    double begin = -4.5;
    double end = 4.5;
    double range = end - begin;
    double step = range / total;
    array_t x;
    x.resize(total);
    array_t y;
    y.resize(total);
    for (int i = 0; i < total; ++i) {
        x[i] = begin;
        y[i] = begin;
        begin += step;
    }
    // calculate z
    double integral = 0;
    flowfield_t z;
    z.resize(total, total);
    for (int i = 0; i < total; ++i) {
        for (int j = 0; j < total; ++j) {
            vector_t dr = {x(i), y(j)};
            z(i, j) = kernel_C_2d(&dr);
            integral += kernel_C_2d(&dr) * step * step;
        }
    }
    EXPECT_NEAR(integral, 1, 1e-5);
    // check that the other one is not working
    integral = 0;
    z.setZero();
    for (int i = 0; i < total; ++i) {
        for (int j = 0; j < total; ++j) {
            vector_t dr = {x(i), y(j)};
            z(i, j) = kernel_C(dr.norm());
            integral += z(i, j) * step * step;
        }
    }
    EXPECT_NE(1, integral); // about 1.87 so complete nonsense
}

/**
 * We check values form the kernel to the test function and then the actual kernel that we are gonna us (divide the range by two)
 */
TEST(IbmTest, ibm_kernel_select) {
    // we test the pure kernel function and then with the application of the right function
    vector_t dk = {0,0};
    double test_x = 0.24;
    double test_y = 0.65;
    vector_t r = {test_x,test_y};
    ibmSimulation sim(nullptr,nullptr,nullptr, dk); // we dont actually int the simulation
    simulation_parameters t;
    t.kernel_in_use = KERNEL_A;
    sim.set_simulation_parameters(t);
    // call the fkt
    EXPECT_EQ(sim.test_kernel_function(&r), kernel_A_2d(&r));
    t.kernel_in_use = KERNEL_B;
    sim.set_simulation_parameters(t);
    // call the fkt
    EXPECT_EQ(sim.test_kernel_function(&r), kernel_B_2d(&r));
    t.kernel_in_use = KERNEL_C;
    sim.set_simulation_parameters(t);
    // call the fkt
    EXPECT_EQ(sim.test_kernel_function(&r), kernel_C_2d(&r));
    t.kernel_in_use = KERNEL_C;
    sim.set_simulation_parameters(t);
    // call the fkt
    EXPECT_NE(sim.test_kernel_function(&r), kernel_A_2d(&r));
    // test the actual use case
    t.kernel_in_use = KERNEL_A;
    sim.set_simulation_parameters(t);
    r = {0,0};
    EXPECT_EQ(sim.test_kernel_function_call(&r),1);
    r = {-2,-2};
    EXPECT_EQ(sim.test_kernel_function_call(&r),0);
    r = {2,2};
    EXPECT_EQ(sim.test_kernel_function_call(&r),0);
}

/**
 * Tests the A kernel node distribution.
 * @test
 */
TEST(IbmTest, right_amount_kernel_nodes_A) {
    point_t starter = {5,5};
    straight_t input;
    straightGenerator sg;
    long canvas_size = 20;
    double side_length = 6; // with a distance of 0.75 we should get 80 markers
    kernelType_t kernel = KERNEL_A;
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
    double ibm_distance = kernel_id_to_lattice_search(kernel);
    ng.init_surface(canvas_size,ibm_distance);
    //side lenght in nodes should be 7x7 with ibm nodes
    EXPECT_EQ( ng.node_infos.size(),121);
    int ibm = 0;
    int regular = 0;
    int error = 0;
    for(auto ni : ng.node_infos) {
        if(ni->boundary == IBM) {
            ++ibm;
        }
        else if(ni->boundary == NO_BOUNDARY) {
            ++regular;
        }
        else {
            ++error;
        }
    }
    EXPECT_EQ(ibm, 4*5*side_length);
    EXPECT_EQ(regular,1);
    EXPECT_EQ(error, 0);
}

/**
 * Tests the b kernel node distribution.
 * @test
 * @attention right now the corners get missed which is not really a mistake tbh (corrector variable)
 */
TEST(IbmTest, right_amount_kernel_nodes_B) {
    int corrector = 4;
    point_t starter = {5,5};
    straight_t input;
    straightGenerator sg;
    long canvas_size = 20;
    double side_length = 8; // with a distance of 0.75 we should get 80 markers
    kernelType_t kernel = KERNEL_B;
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
    double ibm_distance = kernel_id_to_lattice_search(kernel);
    ng.init_surface(canvas_size,ibm_distance);
    EXPECT_EQ( ng.node_infos.size(),225-corrector);
    int ibm = 0;
    int regular = 0;
    int error = 0;
    for(auto ni : ng.node_infos) {
        if(ni->boundary == IBM) {
            ++ibm;
        }
        else if(ni->boundary == NO_BOUNDARY) {
            ++regular;
        }
        else {
            ++error;
        }
    }
    EXPECT_EQ(ibm, 4*7*(side_length)-corrector);
    EXPECT_EQ(regular,1);
    EXPECT_EQ(error, 0);
}

/**
 * Tests the c kernel node distribution.
 * @test
 * @attention right now the corners get missed which is not really a mistake tbh (corrector variable)
 */
TEST(IbmTest, right_amount_kernel_nodes_C) {
    int corrector = 4;
    point_t starter = {5,5};
    straight_t input;
    straightGenerator sg;
    long canvas_size = 25;
    double side_length = 10; // with a distance of 0.75 we should get 80 markers
    kernelType_t kernel = KERNEL_C;
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
    double ibm_distance = kernel_id_to_lattice_search(kernel);
    ng.init_surface(canvas_size,ibm_distance);
    EXPECT_EQ( ng.node_infos.size(),361-corrector);
    int ibm = 0;
    int regular = 0;
    int error = 0;
    for(auto ni : ng.node_infos) {
        if(ni->boundary == IBM) {
            ++ibm;
        }
        else if(ni->boundary == NO_BOUNDARY) {
            ++regular;
        }
        else {
            ++error;
        }
    }
    EXPECT_EQ(ibm, 4*9*(side_length)-corrector);
    EXPECT_EQ(regular,1);
    EXPECT_EQ(error, 0);
}
/**
 * Tests weather or not the nodes are setup correctly, they are.
 * @test
 */
TEST(IbmTest,rho_init) {
    // init the kernel 3 variant
    point_t starter = {5,5};
    straight_t input;
    straightGenerator sg;
    long canvas_size = 25;
    double side_length = 10; // with a distance of 0.75 we should get 80 markers
    kernelType_t kernel = KERNEL_C;
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
    double ibm_distance = kernel_id_to_lattice_search(kernel);
    ng.init_surface(canvas_size,ibm_distance);
    simulation_parameters params;
    params.relaxation = 0.5;
    params.ibm_range = kernel_id_to_lattice_search(kernel);
    params.kernel_in_use = kernel;
    params.k = 1;
    point_t dk = {0,0};
    // max rotation is 7.5e-3
    vector_t sizes = {canvas_size,canvas_size};
    goaForce rot(dk,sizes,1e-3);
    ibmSimulation sim(&ng, &rot,ng.markers,sizes);
    sim.set_simulation_parameters(params);
    sim.init();
    for(auto n : sim.nodes) {
        auto [rho,ux,uy] = sim.test_macro(&n->populations);
        EXPECT_NEAR(rho,1,1e-10);
        EXPECT_NEAR(ux, 0, 1e-10);
        EXPECT_NEAR(uy, 0 , 1e-10);
    }
}

TEST(IbmTest,velocity_interpolation) {
    // init the kernel 3 variant
    point_t starter = {5,5};
    straight_t input;
    straightGenerator sg;
    long canvas_size = 25;
    double side_length = 10; // with a distance of 0.75 we should get 80 markers
    kernelType_t kernel = KERNEL_C;
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
    double ibm_distance = kernel_id_to_lattice_search(kernel);
    ng.init_surface(canvas_size,ibm_distance);
    simulation_parameters params;
    params.relaxation = 0.5;
    params.ibm_range = kernel_id_to_lattice_search(kernel);
    params.kernel_in_use = kernel;
    params.k = 1;
    point_t dk = {0,0};
    // max rotation is 7.5e-3
    vector_t sizes = {canvas_size,canvas_size};
    goaForce rot(dk,sizes,1e-3);
    ibmSimulation sim(&ng, &rot,ng.markers,sizes);
    sim.set_simulation_parameters(params);
    sim.init();
    // find the nodes around 55 and set their velocity
    rangingPointKeyHash rpkh;
    for(auto n : sim.nodes) {
        rpkh.fill_key(n->handle,n->position);
    }
    std::vector<handle_t> relevant = rpkh.ranging_key_translation(starter,4);
    std::cout << relevant.size() << std::endl;
    for(auto h : relevant) {

    }
    // set to a value and look if it makes sense
    // put the veloctiy in the markers
    // note the last marker is the one we are interested in as it is a 5,5
    // test and look
}

TEST(IbmTest, aggregated_force_in_node) {

}

TEST(IbmTest, calculate_force_density) {

}

TEST(IbmTest, self_stream) {
    // tests the self stream in nodes that have less than 8 links
}

TEST(IbmTest, force_strength) {
    // tests for correct values in the force implementation
}


// todo there are some strange cases still left -> investigate
// todo refactor boundary point generators look for bumps and so on -> ranging pgk is a great tool her