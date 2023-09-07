// Include a library file to make sure proper includes are set
#include "simulation.h"
#include "two_step_simulation.h"
#include "one_step_simulation.h"
#include "node.h"
#include "functions.h"
#include "helper_functions.h"
#include "neighborhood.h"
#include "image_converter.h"
#include <gtest/gtest.h>

/**
 *  Tests the streaming in channel 1.
 *  @test
 *  @see basicSimulation::streaming_step_1()
 *  @see basicSimulation::streaming_step_2()
 */
TEST(StreamTests, one_D_streaming_channel_one) {
    // for the streaming test the last two points have to be ignored (aka DRY nodes)
    int size = 20;
    point_t start = {0,0};
    point_t end = {size,0};
    point_t sim_area = {size,1};
    boundaryPointConstructor boundaries(sim_area);
    boundaries.init_structure();
    boundaries.set_point(&start,NO_BOUNDARY);
    boundaries.set_point(&end  ,NO_BOUNDARY);
    basicSimulation sm(&boundaries);
    sm.init();
    // zero the data
    for(auto node : sm.nodes) {
        node->population_even.setZero();
        node->population_odd.setZero();
    }
    //manual streaming channel 1 start at handle 1 at 1,0 to handle 4 at 4,0
    sm.nodes.at(0)->population_even(1) = 1;
    for(int i = 1; i < size; ++i) {
        // go check
        for(auto node : sm.nodes) {
            int expect = 0;
            if(node->handle == i) {
                expect = 1;
            }
            EXPECT_EQ(expect, node->population_even(1));
        }
        // execute two steps
        sm.streaming_step_1();
        sm.streaming_step_2();
    }
}
/**
 *  Tests the streaming in channel 3.
 *  @test
 *  @see basicSimulation::streaming_step_1()
 *  @see basicSimulation::streaming_step_2()
 */
TEST(StreamTests, one_D_streaming_channel_three) {
    int size = 20;
    point_t start = {0,0};
    point_t end = {size,0};
    point_t sim_area = {size,1};
    boundaryPointConstructor boundaries(sim_area);
    boundaries.init_structure();
    boundaries.set_point(&start,NO_BOUNDARY);
    boundaries.set_point(&end  ,NO_BOUNDARY);
    basicSimulation sm(&boundaries);
    sm.init();
    // zero the data
    for(auto node : sm.nodes) {
        node->population_even.setZero();
        node->population_odd.setZero();
    }
    //manual streaming channel 1
    // the -3 is wierd, but if u think about it last tow nodes
    // have handle 20 and 19 in array 19 and 18 so last on is 17 in array
    sm.nodes.at(size-3)->population_even(3) = 1;
    // aka 20 nodes in total numbered 1 to 20 last two a boundary so we ignore those
    for(int i = size-2; i > 0; --i) {
        // go check
        for(auto node : sm.nodes) {
            // std::cout << node->neighbors.size() << std::endl;
            int expect = 0;
            if(node->handle == i) {
                expect = 1;
            }
            EXPECT_EQ(expect, node->population_even(3));
        }
        // execute two steps
        sm.streaming_step_1();
        sm.streaming_step_2();
    }
}

/**
 *  Tests the streaming in channel 2.
 *  @test
 *  @see basicSimulation::streaming_step_1()
 *  @see basicSimulation::streaming_step_2()
 */
TEST(StreamTests, one_D_streaming_channel_two) {
    int size = 20;
    int channel = 2;
    point_t start = {0,0};
    point_t end = {0,size};
    point_t sim_area = {1,size};
    // need to be able to prime the boundary point construtor
    vector_t node_generation = {0,1};
    boundaryPointConstructor boundaries(sim_area);
    boundaries.init_structure();
    boundaries.set_point(&start,NO_BOUNDARY);
    boundaries.set_point(&end  ,NO_BOUNDARY);
    nodeGenerator gen(&boundaries);
    gen.set_discovery_vector(node_generation);
    gen.init();
    basicSimulation sm(&boundaries,&gen);
    sm.init();
    // zero the data
    for(auto node : sm.nodes) {
        node->population_even.setZero();
        node->population_odd.setZero();
    }
    //manual streaming channel 1
    // the -3 is wierd, but if u think about it last tow nodes
    // have handle 20 and 19 in array 19 and 18 so last on is 17 in array
    sm.nodes.at(0)->population_even(channel) = 1;
    // aka 20 nodes in total numbered 1 to 20 last two a boundary so we ignore those
    for(int i = 1; i < size; ++i) {
        // go check
        for(auto node : sm.nodes) {
            // std::cout << node->neighbors.size() << std::endl;
            int expect = 0;
            if(node->handle == i) {
                expect = 1;
            }
            EXPECT_EQ(expect, node->population_even(channel));
        }
        // execute two steps
        sm.streaming_step_1();
        sm.streaming_step_2();
    }
}

/**
 *  Tests the streaming in channel 4.
 *  @test
 *  @see basicSimulation::streaming_step_1()
 *  @see basicSimulation::streaming_step_2()
 */
TEST(StreamTests, one_D_streaming_channel_four) {
    int size = 20;
    int channel = 4;
    point_t start = {0,0};
    point_t end = {0,size};
    point_t sim_area = {1,size};
    // need to be able to prime the boundary point construtor
    vector_t node_generation = {0,1};
    boundaryPointConstructor boundaries(sim_area);
    boundaries.init_structure();
    boundaries.set_point(&start,NO_BOUNDARY);
    boundaries.set_point(&end  ,NO_BOUNDARY);
    nodeGenerator gen(&boundaries);
    gen.set_discovery_vector(node_generation);
    gen.init();
    basicSimulation sm(&boundaries,&gen);
    sm.init();
    // zero the data
    for(auto node : sm.nodes) {
        node->population_even.setZero();
        node->population_odd.setZero();
    }
    //manual streaming channel 1
    // the -3 is wierd, but if u think about it last tow nodes
    // have handle 20 and 19 in array 19 and 18 so last on is 17 in array
    sm.nodes.at(size-3)->population_even(channel) = 1;
    // aka 20 nodes in total numbered 1 to 20 last two a boundary so we ignore those
    for(int i = size-2; i > 0; --i) {
        // go check
        for(auto node : sm.nodes) {
            // std::cout << node->neighbors.size() << std::endl;
            int expect = 0;
            if(node->handle == i) {
                expect = 1;
            }
            EXPECT_EQ(expect, node->population_even(channel));
        }
        // execute two steps
        sm.streaming_step_1();
        sm.streaming_step_2();
    }
}

/**
 *  Tests the streaming in channel 5.
 *  @test
 *  @see basicSimulation::streaming_step_1()
 *  @see basicSimulation::streaming_step_2()
 */
TEST(StreamTests, one_D_streaming_channel_five) {
    int size = 20;
    int channel = 5;
    point_t start = {0,0};
    point_t end = {size,size};
    point_t sim_area = {size,size};
    // need to be able to prime the boundary point construtor
    vector_t node_generation = {1,1};
    boundaryPointConstructor boundaries(sim_area);
    boundaries.init_structure();
    boundaries.set_point(&start,NO_BOUNDARY);
    boundaries.set_point(&end  ,NO_BOUNDARY);
    nodeGenerator gen(&boundaries);
    gen.set_discovery_vector(node_generation);
    gen.init();
    basicSimulation sm(&boundaries,&gen);
    sm.init();
    // zero the data
    for(auto node : sm.nodes) {
        node->population_even.setZero();
        node->population_odd.setZero();
    }
    //manual streaming channel 1
    // the -3 is wierd, but if u think about it last tow nodes
    // have handle 20 and 19 in array 19 and 18 so last on is 17 in array
    sm.nodes.at(0)->population_even(channel) = 1;
    // aka 20 nodes in total numbered 1 to 20 last two a boundary so we ignore those
    for(int i = 1; i < size; ++i) {
        // go check
        for(auto node : sm.nodes) {
            // std::cout << node->neighbors.size() << std::endl;
            int expect = 0;
            if(node->handle == i) {
                expect = 1;
            }
            EXPECT_EQ(expect, node->population_even(channel));
        }
        // execute two steps
        sm.streaming_step_1();
        sm.streaming_step_2();
    }
}

/**
 *  Tests the streaming in channel 7.
 *  @test
 *  @see basicSimulation::streaming_step_1()
 *  @see basicSimulation::streaming_step_2()
 */
TEST(StreamTests, one_D_streaming_channel_seven) {
    int size = 20;
    int channel = 7;
    point_t start = {0,0};
    point_t end = {size,size};
    point_t sim_area = {size,size};
    // need to be able to prime the boundary point construtor
    vector_t node_generation = {1,1};
    boundaryPointConstructor boundaries(sim_area);
    boundaries.init_structure();
    boundaries.set_point(&start,NO_BOUNDARY);
    boundaries.set_point(&end  ,NO_BOUNDARY);
    nodeGenerator gen(&boundaries);
    gen.set_discovery_vector(node_generation);
    gen.init();
    basicSimulation sm(&boundaries,&gen);
    sm.init();
    // zero the data
    for(auto node : sm.nodes) {
        node->population_even.setZero();
        node->population_odd.setZero();
    }
    //manual streaming channel 1
    // the -3 is wierd, but if u think about it last tow nodes
    // have handle 20 and 19 in array 19 and 18 so last on is 17 in array
    sm.nodes.at(size-3)->population_even(channel) = 1;
    // aka 20 nodes in total numbered 1 to 20 last two a boundary so we ignore those
    for(int i = size-2; i > 0; --i) {
        // go check
        for(auto node : sm.nodes) {
            // std::cout << node->neighbors.size() << std::endl;
            int expect = 0;
            if(node->handle == i) {
                expect = 1;
            }
            EXPECT_EQ(expect, node->population_even(channel));
        }
        // execute two steps
        sm.streaming_step_1();
        sm.streaming_step_2();
    }
}

/**
 *  Tests the streaming in channel 6.
 *  @test
 *  @see basicSimulation::streaming_step_1()
 *  @see basicSimulation::streaming_step_2()
 */
TEST(StreamTests, one_D_streaming_channel_six) {
    int size = 20;
    int channel = 6;
    point_t start = {size,0};
    point_t end = {0,size};
    point_t sim_area = {size,size};
    // need to be able to prime the boundary point construtor
    vector_t node_generation = {-1,1};
    boundaryPointConstructor boundaries(sim_area);
    boundaries.init_structure();
    boundaries.set_point(&start,NO_BOUNDARY);
    boundaries.set_point(&end  ,NO_BOUNDARY);
    nodeGenerator gen(&boundaries);
    gen.set_discovery_vector(node_generation);
    gen.init();
    basicSimulation sm(&boundaries,&gen);
    sm.init();
    // zero the data
    for(auto node : sm.nodes) {
        node->population_even.setZero();
        node->population_odd.setZero();
    }
    //manual streaming channel 1
    // the -3 is wierd, but if u think about it last tow nodes
    // have handle 20 and 19 in array 19 and 18 so last on is 17 in array
    sm.nodes.at(0)->population_even(channel) = 1;
    // aka 20 nodes in total numbered 1 to 20 last two a boundary so we ignore those
    for(int i = 1; i < size; ++i) {
        // go check
        for(auto node : sm.nodes) {
            // std::cout << node->neighbors.size() << std::endl;
            int expect = 0;
            if(node->handle == i) {
                expect = 1;
            }
            EXPECT_EQ(expect, node->population_even(channel));
        }
        // execute two steps
        sm.streaming_step_1();
        sm.streaming_step_2();
    }
}

/**
 *  Tests the streaming in channel 8.
 *  @test
 *  @see basicSimulation::streaming_step_1()
 *  @see basicSimulation::streaming_step_2()
 */
TEST(StreamTests, one_D_streaming_channel_eight) {
    int size = 20;
    int channel = 8;
    point_t start = {size,0};
    point_t end = {0,size};
    point_t sim_area = {size,size};
    // need to be able to prime the boundary point construtor
    vector_t node_generation = {-1,1};
    boundaryPointConstructor boundaries(sim_area);
    boundaries.init_structure();
    boundaries.set_point(&start,NO_BOUNDARY);
    boundaries.set_point(&end  ,NO_BOUNDARY);
    nodeGenerator gen(&boundaries);
    gen.set_discovery_vector(node_generation);
    gen.init();
    basicSimulation sm(&boundaries,&gen);
    sm.init();
    // zero the data
    for(auto node : sm.nodes) {
        node->population_even.setZero();
        node->population_odd.setZero();
    }
    //manual streaming channel 1
    // the -3 is wierd, but if u think about it last tow nodes
    // have handle 20 and 19 in array 19 and 18 so last on is 17 in array
    sm.nodes.at(size-3)->population_even(channel) = 1;
    // aka 20 nodes in total numbered 1 to 20 last two a boundary so we ignore those
    for(int i = size-2; i > 0; --i) {
        // go check
        for(auto node : sm.nodes) {
            // std::cout << node->neighbors.size() << std::endl;
            int expect = 0;
            if(node->handle == i) {
                expect = 1;
            }
            EXPECT_EQ(expect, node->population_even(channel));
        }
        // execute two steps
        sm.streaming_step_1();
        sm.streaming_step_2();
    }
}

/**
 *  Tests with specific codes weather or not streaming and collision is consistent.
 *  @test
 *  @see basicSimulation::run()
 */
TEST(StreamTests, combinded_test_boundary_consistent) {
    // checks leakages form the boundary nodes
    // setup sim
    int size = 5;
    int sim_size = size-2;
    int steps = 10;
    point_t p = {size,size};
    boundaryPointConstructor boundaries(p);
    boundaries.init_quader();
    basicSimulation sim(&boundaries);
    sim.init();
    // zero the sim space and
    // set specific codes to the corners and observe them
    // (if someone gets where they from expect spaghetti lmao)
    std::vector<array_t> index_corners;
    std::vector<int> index_codes;
    array_t  corner_1(2); corner_1 << 1,1;
    int code_104 = 104;
    array_t  corner_2(2); corner_2 << 3,1; ;
    int code_107 = 107;
    array_t  corner_3(2); corner_3 << 3,3;
    int code_126 = 126;
    array_t  corner_4(2); corner_4 << 1,3;
    int code_202 = 202;
    index_corners.push_back(corner_1);
    index_corners.push_back(corner_2);
    index_corners.push_back(corner_3);
    index_corners.push_back(corner_4);
    index_codes.push_back(code_104);
    index_codes.push_back(code_107);
    index_codes.push_back(code_126);
    index_codes.push_back(code_202);
    // zero the nodes and then put the codes in
    for(auto node : sim.nodes) {
        if(node->node_type == WET) {
            node->population_even= 0;
            node->population_odd = 0;
            macro(node);
        }
        else {
            node->population_even= 1;
            node->population_odd = 1;
            macro(node);
        }
    }
    // 4 corners
    for( int i = 0; i < 4; ++i) {
        for(auto node: sim.nodes) {
            if(node_position_comparison(node,&index_corners.at(i))) {
                // increment
                for(int j = 0; j < node->population_even.size(); ++j) {
                    node->population_even(j) = index_codes.at(i) * 10 + j;
                }
                node->population_odd = node->population_even;
            }
        }
    }
    // run checks rest is preamble
    for(int i = 0; i < steps; i++) {
        sim.run();
        for(auto node: sim.nodes) {
            if(node->node_type == DRY) {
                EXPECT_EQ(node->rho, 9);
            }
        }
    }
}

/**
 * Tests weather or not bounce back works form channel one to three.
 * @test
 * @see basicSimulation::bounce_back()
 * @see basicSimulation::streaming_step_1()
 * @see basicSimulation::streaming_step_2()
 */
TEST(BounceBackTesting, Horizontals_one_three) {
    // we need to generate two nodes and only really consider the bb
    // generate two points we have two bounary nodes and one wet node
    // in this case we are interested in channel 1 to 3 bb
    int size = 3;
    point_t start = {0,0};
    point_t end = {size,0};
    point_t sim_area = {size,1};
    boundaryPointConstructor boundaries(sim_area);
    boundaries.init_structure();
    boundaries.set_point(&start,BOUNCE_BACK);
    boundaries.set_point(&end  ,BOUNCE_BACK);
    nodeGenerator gen(&boundaries);
    gen.init();
    basicSimulation sm(&boundaries,&gen);
    sm.init();
    // zero the data
    for(auto node : sm.nodes) {
        node->population_even.setZero();
        node->population_odd.setZero();
        /* useful debug fragement
        std::cout << node->handle << std::endl;
        std::cout << node->neighbors.size() << std::endl;
        std::cout << node->position << std::endl;
        std::cout << std::endl;
         */
    }
    // set the value of channel 1 in the middle node to 1
    // should reappear in channel 3 and vice versa
    // we use two points with handle 1 & 2 to check the behaviour of bb in channel 1 and 3
    sm.nodes.at(0)->population_even(3) = 1;
    sm.nodes.at(1)->population_even(1) = 1;
    sm.streaming_step_1();
    sm.streaming_step_2();
    sm.bounce_back();
    EXPECT_EQ(sm.nodes.at(1)->population_even(3),1);
    EXPECT_EQ(sm.nodes.at(0)->population_even(1),1);
}

/**
 * Tests weather or not bounce back works form channel two to four.
 * @test
 * @see basicSimulation::bounce_back()
 * @see basicSimulation::streaming_step_1()
 * @see basicSimulation::streaming_step_2()
 */
TEST(BounceBackTesting, Horizontals_two_four) {
    // we need to generate two nodes and only really consider the bb
    // generate two points we have two bounary nodes and one wet node
    // in this case we are interested in channel 1 to 3 bb
    // need to be able to prime the boundary point construtor
    int size = 3;
    point_t start = {0,0};
    point_t end = {size,size};
    point_t sim_area = {size,size};
    vector_t node_generation = {0,1};
    boundaryPointConstructor boundaries(sim_area);
    boundaries.init_structure();
    boundaries.set_point(&start,BOUNCE_BACK);
    boundaries.set_point(&end  ,BOUNCE_BACK);
    nodeGenerator gen(&boundaries);
    gen.set_discovery_vector(node_generation);
    gen.init();
    basicSimulation sm(&boundaries,&gen);
    sm.init();
    // zero the data
    for(auto node : sm.nodes) {
        node->population_even.setZero();
        node->population_odd.setZero();
        /* useful debug fragement
        std::cout << node->handle << std::endl;
        std::cout << node->neighbors.size() << std::endl;
        std::cout << node->position << std::endl;
        std::cout << std::endl;
         */
    }
    // set the value of channel 1 in the middle node to 1
    // should reappear in channel 3 and vice versa
    // we use two points with handle 1 & 2 to check the behaviour of bb in channel 1 and 3
    sm.nodes.at(0)->population_even(2) = 1;
    sm.nodes.at(1)->population_even(4) = 1;
    sm.streaming_step_1();
    sm.streaming_step_2();
    sm.bounce_back();
    EXPECT_EQ(sm.nodes.at(1)->population_even(2),1);
    EXPECT_EQ(sm.nodes.at(0)->population_even(4),1);
}

/**
 * Tests weather or not bounce back works form channel five to seven.
 * @test
 * @see basicSimulation::bounce_back()
 * @see basicSimulation::streaming_step_1()
 * @see basicSimulation::streaming_step_2()
 */
TEST(BounceBackTesting, Oblique_five_seven) {
    // we need to generate two nodes and only really consider the bb
    // generate two points we have two bounary nodes and one wet node
    // in this case we are interested in channel 1 to 3 bb
    int size = 3;
    point_t start = {0,0};
    point_t end = {size,size};
    point_t sim_area = {size,size};
    // need to be able to prime the boundary point construtor
    vector_t node_generation = {1,1};
    boundaryPointConstructor boundaries(sim_area);
    boundaries.init_structure();
    boundaries.set_point(&start,BOUNCE_BACK);
    boundaries.set_point(&end  ,BOUNCE_BACK);
    nodeGenerator gen(&boundaries);
    gen.set_discovery_vector(node_generation);
    gen.init();
    basicSimulation sm(&boundaries,&gen);
    sm.init();
    // zero the data
    for(auto node : sm.nodes) {
        node->population_even.setZero();
        node->population_odd.setZero();
        //useful debug fragement
        /*
        std::cout << node->handle << std::endl;
        std::cout << node->neighbors.size() << std::endl;
        std::cout << node->position << std::endl;
        std::cout << std::endl;
         */
    }
    // set the value of channel 1 in the middle node to 1
    // should reappear in channel 3 and vice versa
    // we use two points with handle 1 & 2 to check the behaviour of bb in channel 1 and 3
    sm.nodes.at(0)->population_even(7) = 1;
    sm.nodes.at(1)->population_even(5) = 1;
    sm.streaming_step_1();
    sm.streaming_step_2();
    sm.bounce_back();
    EXPECT_EQ(sm.nodes.at(0)->population_even(5),1);
    EXPECT_EQ(sm.nodes.at(1)->population_even(7),1);
}

/**
 * Tests weather or not bounce back works form channel six to eight.
 * @test
 * @see basicSimulation::bounce_back()
 * @see basicSimulation::streaming_step_1()
 * @see basicSimulation::streaming_step_2()
 */
TEST(BounceBackTesting, Oblique_six_eight) {
    // we need to generate two nodes and only really consider the bb
    // generate two points we have two bounary nodes and one wet node
    // in this case we are interested in channel 1 to 3 bb
    int size = 3;
    point_t start = {0, size};
    point_t end = {size, 0};
    point_t sim_area = {size, size};
    // need to be able to prime the boundary point construtor
    vector_t node_generation = {1, -1};
    boundaryPointConstructor boundaries(sim_area);
    boundaries.init_structure();
    boundaries.set_point(&start, BOUNCE_BACK);
    boundaries.set_point(&end, BOUNCE_BACK);
    nodeGenerator gen(&boundaries);
    gen.set_discovery_vector(node_generation);
    gen.init();
    basicSimulation sm(&boundaries, &gen);
    sm.init();
    // zero the data
    for (auto node : sm.nodes) {
        node->population_even.setZero();
        node->population_odd.setZero();
    }
    // set the value of channel 1 in the middle node to 1
    // should reappear in channel 3 and vice versa
    // we use two points with handle 1 & 2 to check the behaviour of bb in channel 1 and 3
    sm.nodes.at(0)->population_even(6) = 1;
    sm.nodes.at(1)->population_even(8) = 1;
    sm.streaming_step_1();
    sm.streaming_step_2();
    sm.bounce_back();
    EXPECT_EQ(sm.nodes.at(0)->population_even(8), 1);
    EXPECT_EQ(sm.nodes.at(1)->population_even(6), 1);
}

/**
 * Tests out weather or not a moving bounce back boundary has the right result.
 * @test
 * @see basicSimulation::bounce_back()
 * @see basicSimulation::streaming_step_1()
 * @see basicSimulation::streaming_step_2()
 */
TEST(BounceBackTesting, moving) {
    // test the actual moving code
    // if this still doesnt work out ill do poisioulle flow...
    // probable cause leakage into channels that have nothing to do
    // in the corners ?!
    // put a bunch of 1 into the top middle node of a 3x3 simspace
    int size = 5;
    int sim_size = size-2;
    double uw = 0.1;
    point_t p = {size,size};
    boundaryPointConstructor boundaries(p);
    boundaries.init_sliding_lid();
    basicSimulation sim(&boundaries);
    sim.init();
    double re = 1000;
    double base_length = sim_size - 2;
    simulation_parameters params;
    params.u_wall = 0.1;
    params.relaxation = (2*re)/(6*base_length*params.u_wall+re);
    sim.set_simulation_parameters(params);
    // zero the data
    for(auto node : sim.nodes) {
        node->population_even.setZero();
        node->population_odd.setZero();
    }
    // run the bb sequence
    sim.streaming_step_1();
    sim.streaming_step_2();
    sim.bounce_back();
    // check results
    double channel7 = -1.0/6* uw;
    double channel8 = 1.0/6 *uw;
    for(auto node : sim.nodes) {
        if(node->node_type == WET) {
            if(node->position(1) == sim_size) {
                EXPECT_EQ(node->population_even(7),channel7);
                EXPECT_EQ(node->population_even(8),channel8);
            }
        }
    }
}

/**
 * Tests if the values in channel 0 reamain untouched by streaming.
 * @test
 * @see basicSimulation::streaming_step_1()
 * @see basicSimulation::streaming_step_2()
 */
TEST(StreamTests, channel_0_persistent) {
    // test weather or not values in channel 0 persist throu streaming
    int size = 5;
    int steps = 5;
    point_t p = {size,size};
    boundaryPointConstructor boundaries(p);
    boundaries.init_sliding_lid();
    // simulation init
    basicSimulation sim(&boundaries);
    sim.init();
    // parameters
    double re = 1000;
    double base_length = size - 2;
    simulation_parameters params;
    params.u_wall = 0.1;
    params.relaxation = (2*re)/(6*base_length*params.u_wall+re);
    sim.set_simulation_parameters(params);

    //
    for(int i = 0; i < steps; ++i) {
        double before_data_value = sim.nodes.at(0)->population_even(0);
        sim.streaming_step_1();
        sim.streaming_step_2();
        sim.bounce_back();
        // after streaming the before data value in channel 0 should still be there
        EXPECT_EQ(before_data_value, sim.nodes.at(0)->population_even(0));
        for(auto n : sim.nodes) {
            macro(n);
        }
        for(auto n : sim.nodes) {
            collision(n,params.relaxation);
        }

    }
}

/**
 * Generates 2 node generators one that writes to file and the other one not.
 * Checks weather or not results are the same.
 * Tests the write and readback.
 * @test
 * @see nodeGenerator::set_redo_save()
 */
TEST(FunctionalTest, read_write) {
    int size = 5;
    point_t p = {size,size};
    boundaryPointConstructor boundaries(p);
    boundaries.init_quader();
    // run the first node generator
    nodeGenerator nodes(&boundaries);
    nodes.set_redo_save(true, true);
    nodes.init();
    // run a second one but relay on the result of the first one
    nodeGenerator nodes2(&boundaries);
    nodes2.set_redo_save(false,false);
    nodes2.init();
    // check for right amout of nodes in both containers
    EXPECT_EQ(nodes.node_infos.size(),size*size );
    ASSERT_EQ(nodes.node_infos.size(),nodes2.node_infos.size());
    // check the values they should agree
    for(int i = 0; i < nodes.node_infos.size(); ++i) {
        auto n1 = nodes.node_infos.at(i);
        auto n2 = nodes2.node_infos.at(i);
        // check
        EXPECT_EQ(n1->handle,n2->handle);
        EXPECT_EQ(n1->type,n2->type);
        EXPECT_EQ(n1->boundary,n2->boundary);
        // check prior sizes
        ASSERT_EQ(n1->position.size(), n2->position.size());
        ASSERT_EQ(n1->links.size(), n2->links.size());
        for(int j = 0; j < n1->position.size(); ++j) {
            EXPECT_EQ(n1->position(j), n2->position(j));
        }
        for(int j = 0; j < n1->links.size(); ++j) {
            EXPECT_EQ(n1->links.at(j).handle, n2->links.at(j).handle);
            EXPECT_EQ(n1->links.at(j).channel,n2->links.at(j).channel);
        }
    }
}

/**
 * Tests the functionaltiy of one step streaming in channel 1 and 3.
 * @test
 * @see basicSimulation::fused_streaming()
 */
TEST(StreamTests, fused_streaming_13) {
    // sizes matters
    int size = 4;
    point_t start = {0,0};
    point_t end = {size,0};
    point_t sim_area = {size,size-1};
    // init
    boundaryPointConstructor boundaries(sim_area);
    boundaries.init_quader();
    // boundaries.visualize_2D_boundary();
    nodeGenerator gen(&boundaries);
    gen.init_fused(size);
    // gen.visualize_2D_nodes();
    EXPECT_EQ(gen.node_infos.size(),2);
    // gen.visualize_2D_nodes(4);
    basicSimulation sm(&boundaries, &gen);
    sm.fused_init();
    // check sim sizes
    EXPECT_EQ(sm.nodes.size(),2);
    // zero data
    // zero the data
    for(auto node : sm.nodes) {
        node->population_even.setZero();
        node->population_odd.setZero();
    }
    // set the value of channel 1 in the middle node to 1
    // should reappear in channel 3 and vice versa
    // we use two points with handle 1 & 2 to check the behaviour of bb in channel 1 and 3
    sm.nodes.at(0)->population_even(3) = 1;
    sm.nodes.at(1)->population_even(1) = 1;
    // do steps and check correct positions
    for(auto n : sm.nodes) {
        sm.fused_streaming(n);
    }
    for(auto n : sm.nodes) {
        // switchero
        array_t * temp = n->current_population;
        n->current_population = n->next_population;
        n->next_population = temp;
    }
    // channels switched (bb)
    EXPECT_EQ(sm.nodes.at(0)->current_population->operator()(1), 1);
    EXPECT_EQ(sm.nodes.at(1)->current_population->operator()(3), 1);
    // propagation test
    for(auto n : sm.nodes) {
        sm.fused_streaming(n);
    }
    for(auto n : sm.nodes) {
        // switchero
        array_t * temp = n->current_population;
        n->current_population = n->next_population;
        n->next_population = temp;
    }
    EXPECT_EQ(sm.nodes.at(0)->current_population->operator()(3), 1);
    EXPECT_EQ(sm.nodes.at(1)->current_population->operator()(1), 1);
}

/**
 * Tests the functionality of one step streaming in channel 2 and 4.
 * @test
 * @see basicSimulation::fused_streaming()
 */
TEST(StreamTests, fused_streaming_24) {
    // sizes matters
    int size = 4;
    point_t start = {0,0};
    point_t end = {size,0};
    point_t sim_area = {size-1,size};
    // init
    boundaryPointConstructor boundaries(sim_area);
    boundaries.init_quader();
    // boundaries.visualize_2D_boundary();
    nodeGenerator gen(&boundaries);
    gen.init_fused(size);
    // gen.visualize_2D_nodes();
    basicSimulation sm(&boundaries, &gen);
    sm.fused_init();
    // check sim sizes
    EXPECT_EQ(sm.nodes.size(),2);
    // zero data
    // zero the data
    for(auto node : sm.nodes) {
        node->population_even.setZero();
        node->population_odd.setZero();
    }
    // set the value of channel 1 in the middle node to 1
    // should reappear in channel 3 and vice versa
    // we use two points with handle 1 & 2 to check the behaviour of bb in channel 1 and 3
    sm.nodes.at(0)->population_even(4) = 1;
    sm.nodes.at(1)->population_even(2) = 1;
    // do steps and check correct positions
    for(auto n : sm.nodes) {
        sm.fused_streaming(n);
    }
    for(auto n : sm.nodes) {
        // switchero
        array_t * temp = n->current_population;
        n->current_population = n->next_population;
        n->next_population = temp;
    }
    // channels switched (bb)
    EXPECT_EQ(sm.nodes.at(0)->current_population->operator()(2), 1);
    EXPECT_EQ(sm.nodes.at(1)->current_population->operator()(4), 1);
    // propagation test
    for(auto n : sm.nodes) {
        sm.fused_streaming(n);
    }
    for(auto n : sm.nodes) {
        // switchero
        array_t * temp = n->current_population;
        n->current_population = n->next_population;
        n->next_population = temp;
    }
    EXPECT_EQ(sm.nodes.at(0)->current_population->operator()(4), 1);
    EXPECT_EQ(sm.nodes.at(1)->current_population->operator()(2), 1);
}

/**
 * Tests the different implementations of macro and collision against each other.
 * @test
 * @see macro()
 * @see fused_macro()
 * @see collision()
 * @see fused_collision()
 * @see optimizedSimulation::one_step_macro_collision()
 */
TEST(FunctionalTest, one_step_macro_collison) {
    // test against original
    double relaxation_time = 0.5;
    optimizedSimulation oLuv(nullptr, nullptr);
    simulation_parameters simulation_params;
    simulation_params.relaxation = relaxation_time; // tau
    oLuv.set_simulation_parameters(simulation_params);
    // init stuff
    point_t pos = {0,0};
    node node_original(1,velocity_set.rows(),velocity_set.cols(),pos,NO_BOUNDARY);
    node node_fused(1,velocity_set.rows(),velocity_set.cols(),pos,NO_BOUNDARY);
    oNode s_node(1,velocity_set.cols(),NO_BOUNDARY);
    // set good values
    node_original.population_even.setOnes();
    node_fused.population_even.setOnes();
    s_node.populations.setOnes();
    // run functions
    macro(&node_original);
    fused_macro(&node_fused);
    collision(&node_original,relaxation_time);
    fused_collision(&node_fused, relaxation_time);
    oLuv.one_step_macro_collision(&s_node);
    // compare values
    int o= 0;
    for(int i = 0; i < velocity_set.cols(); ++i) {
        EXPECT_EQ(node_original.population_even(i),node_fused.population_even(i));
        EXPECT_EQ(s_node.populations(o + i),node_original.population_even(i));
    }
}

/**
 * Tests the optimized implementation of the streaming in channel 1 and 3.
 * @test
 * @see optimizedSimulation::streaming()
 */
TEST(StreamTests, oSimu_streaming_13) {
    // sizes matters
    int step = 0;
    int size = 4;
    point_t start = {0,0};
    point_t end = {size,0};
    point_t sim_area = {size,size-1};
    // init
    boundaryPointConstructor boundaries(sim_area);
    boundaries.init_quader();
    //boundaries.visualize_2D_boundary(size);
    nodeGenerator gen(&boundaries);
    gen.init_fused(size);
    //gen.visualize_2D_nodes(4);
    optimizedSimulation sm(&boundaries, &gen);
    sm.init();
    // check sim sizes
    EXPECT_EQ(sm.nodes.size(),2);
    // zero data
    // zero the data
    for(auto node : sm.nodes) {
        node->populations.setZero();
    }
    // set the value of channel 1 in the middle node to 1
    // should reappear in channel 3 and vice versa
    // we use two points with handle 1 & 2 to check the behaviour of bb in channel 1 and 3
    sm.nodes.at(0)->populations(3) = 1;
    sm.nodes.at(1)->populations(1) = 1;
    // do steps and check correct positions
    // channels switched (bb)
    sm.offset_sim = ((step +1) & 0x1) * 9;
    sm.offset_node = (step & 0x1) * 9;
    for(auto n : sm.nodes) {
        sm.streaming(n);
    }
    step++;
    EXPECT_EQ(sm.nodes.at(0)->populations(1 + sm.offset_sim), 1);
    EXPECT_EQ(sm.nodes.at(1)->populations(3 + sm.offset_sim), 1);
    // propagation test
    sm.offset_sim = ((step +1) & 0x1) * 9;
    sm.offset_node = (step & 0x1) * 9;
    for(auto n : sm.nodes) {
        sm.streaming(n);
    }
    step++;
    EXPECT_EQ(sm.nodes.at(0)->populations(3 + sm.offset_sim), 1);
    EXPECT_EQ(sm.nodes.at(1)->populations(1 + sm.offset_sim), 1);
}

/**
 * Tests the optimized implementation of the streaming in channel 2 and 4.
 * @test
 * @see optimizedSimulation::streaming()
 */
TEST(StreamTests, oSimu_streaming_24) {
    // sizes matters
    int step = 0;
    int size = 4;
    point_t start = {0,0};
    point_t end = {size,0};
    point_t sim_area = {size - 1,size};
    // init
    boundaryPointConstructor boundaries(sim_area);
    boundaries.init_quader();
    //boundaries.visualize_2D_boundary(size);
    nodeGenerator gen(&boundaries);
    gen.init_fused(size);
    //gen.visualize_2D_nodes(4);
    optimizedSimulation sm(&boundaries, &gen);
    sm.init();
    // check sim sizes
    EXPECT_EQ(sm.nodes.size(),2);
    // zero data
    // zero the data
    for(auto node : sm.nodes) {
        node->populations.setZero();
    }
    // set the value of channel 1 in the middle node to 1
    // should reappear in channel 3 and vice versa
    // we use two points with handle 1 & 2 to check the behaviour of bb in channel 1 and 3
    sm.nodes.at(0)->populations(4) = 1;
    sm.nodes.at(1)->populations(2) = 1;
    // do steps and check correct positions
    // channels switched (bb)
    sm.offset_sim = ((step +1) & 0x1) * 9;
    sm.offset_node = (step & 0x1) * 9;
    for(auto n : sm.nodes) {
        sm.streaming(n);
    }
    step++;
    EXPECT_EQ(sm.nodes.at(0)->populations(2 + sm.offset_sim), 1);
    EXPECT_EQ(sm.nodes.at(1)->populations(4 + sm.offset_sim), 1);
    // propagation test
    sm.offset_sim = ((step +1) & 0x1) * 9;
    sm.offset_node = (step & 0x1) * 9;
    for(auto n : sm.nodes) {
        sm.streaming(n);
    }
    step++;
    EXPECT_EQ(sm.nodes.at(0)->populations(4 + sm.offset_sim), 1);
    EXPECT_EQ(sm.nodes.at(1)->populations(2 + sm.offset_sim), 1);
}

/**
 * Tests the optimized implementation of the streaming in channel 5 and 7.
 * @test
 * @see optimizedSimulation::streaming()
 */
TEST(StreamTests, oSimu_streaming_57) {
    // sizes matters
    int step = 0;
    int size = 4;
    point_t start = {0,0};
    point_t end = {size,0};
    point_t sim_area = {size,size};
    // init boundaries
    boundaryPointConstructor boundaries(sim_area);
    boundaries.init_structure();
    point_t point = {0,0};
    boundaries.set_point(&point,BOUNCE_BACK);
    point = {1,0};
    boundaries.set_point(&point,BOUNCE_BACK);
    point = {2,0};
    boundaries.set_point(&point,BOUNCE_BACK);
    point = {2,1};
    boundaries.set_point(&point,BOUNCE_BACK);
    point = {3,1};
    boundaries.set_point(&point,BOUNCE_BACK);
    point = {3,2};
    boundaries.set_point(&point,BOUNCE_BACK);
    point = {3,3};
    boundaries.set_point(&point,BOUNCE_BACK);
    point = {2,3};
    boundaries.set_point(&point,BOUNCE_BACK);
    point = {1,3};
    boundaries.set_point(&point,BOUNCE_BACK);
    point = {1,2};
    boundaries.set_point(&point,BOUNCE_BACK);
    point = {0,2};
    boundaries.set_point(&point,BOUNCE_BACK);
    point = {0,1};
    boundaries.set_point(&point,BOUNCE_BACK);
    // boundaries.visualize_2D_boundary();
    // node init
    nodeGenerator gen(&boundaries);
    gen.init_fused(size);
    // gen.visualize_2D_nodes(size);
    // tests
    optimizedSimulation sm(&boundaries, &gen);
    sm.init();
    // check sim sizes
    EXPECT_EQ(sm.nodes.size(),2);
    // zero data
    // zero the data
    for(auto node : sm.nodes) {
        node->populations.setZero();
    }
    sm.nodes.at(0)->populations(7) = 1;
    sm.nodes.at(1)->populations(5) = 1;
    // do steps and check correct positions
    // channels switched (bb)
    sm.offset_sim = ((step +1) & 0x1) * 9;
    sm.offset_node = (step & 0x1) * 9;
    for(auto n : sm.nodes) {
        sm.streaming(n);
    }
    step++;
    EXPECT_EQ(sm.nodes.at(0)->populations(5+ sm.offset_sim), 1);
    EXPECT_EQ(sm.nodes.at(1)->populations(7 + sm.offset_sim), 1);
    // propagation test
    sm.offset_sim = ((step +1) & 0x1) * 9;
    sm.offset_node = (step & 0x1) * 9;
    for(auto n : sm.nodes) {
        sm.streaming(n);
    }
    step++;
    EXPECT_EQ(sm.nodes.at(0)->populations(7 + sm.offset_sim), 1);
    EXPECT_EQ(sm.nodes.at(1)->populations(5 + sm.offset_sim), 1);
}

/**
 * Tests the optimized implementation of the streaming in channel 6 and 8.
 * @test
 * @see optimizedSimulation::streaming()
 */
TEST(StreamTests, oSimu_streaming_68) {
    // sizes matters
    int step = 0;
    int size = 4;
    point_t start = {0,0};
    point_t end = {size,0};
    point_t sim_area = {size,size};
    // init boundaries
    boundaryPointConstructor boundaries(sim_area);
    // manual setup
    boundaries.init_structure();
    point_t point = {1,0};
    boundaries.set_point(&point,BOUNCE_BACK);
    point = {2,0};
    boundaries.set_point(&point,BOUNCE_BACK);
    point = {3,0};
    boundaries.set_point(&point,BOUNCE_BACK);
    point = {3,1};
    boundaries.set_point(&point,BOUNCE_BACK);
    point = {3,2};
    boundaries.set_point(&point,BOUNCE_BACK);
    point = {2,2};
    boundaries.set_point(&point,BOUNCE_BACK);
    point = {2,3};
    boundaries.set_point(&point,BOUNCE_BACK);
    point = {1,3};
    boundaries.set_point(&point,BOUNCE_BACK);
    point = {0,3};
    boundaries.set_point(&point,BOUNCE_BACK);
    point = {0,2};
    boundaries.set_point(&point,BOUNCE_BACK);
    point = {0,1};
    boundaries.set_point(&point,BOUNCE_BACK);
    point = {1,1};
    boundaries.set_point(&point,BOUNCE_BACK);
    // need to reorder nodes!!!!!!!! so that the surface is closed!!
    // new init does this
    // boundaries.visualize_2D_boundary();
    // node init
    // first found points seems to be bugged 2.5,2.5
    nodeGenerator gen(&boundaries);
    gen.init_fused(size);
    // gen.visualize_2D_nodes(size);
    optimizedSimulation sm(&boundaries, &gen);
    sm.init();
    // check sim sizes
    EXPECT_EQ(sm.nodes.size(),2);
    // zero data
    // zero the data
    for(auto node : sm.nodes) {
        node->populations.setZero();
    }
    sm.nodes.at(0)->populations(8) = 1;
    sm.nodes.at(1)->populations(6) = 1;
    // do steps and check correct positions
    // channels switched (bb)
    sm.offset_sim = ((step +1) & 0x1) * 9;
    sm.offset_node = (step & 0x1) * 9;
    for(auto n : sm.nodes) {
        sm.streaming(n);
    }
    step++;
    EXPECT_EQ(sm.nodes.at(0)->populations(6+ sm.offset_sim), 1);
    EXPECT_EQ(sm.nodes.at(1)->populations(8 + sm.offset_sim), 1);
    // propagation test
    sm.offset_sim = ((step +1) & 0x1) * 9;
    sm.offset_node = (step & 0x1) * 9;
    for(auto n : sm.nodes) {
        sm.streaming(n);
    }
    step++;
    EXPECT_EQ(sm.nodes.at(0)->populations(8 + sm.offset_sim), 1);
    EXPECT_EQ(sm.nodes.at(1)->populations(6 + sm.offset_sim), 1);
}
