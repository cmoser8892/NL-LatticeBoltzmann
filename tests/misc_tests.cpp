// Include a library file to make sure proper includes are set
#include "simulation.h"
#include "node.h"
#include "functions.h"
#include "helper_functions.h"
#include "neighborhood.h"
#include "image_converter.h"
#include <gtest/gtest.h>

TEST(MiscTest, TestingWhatever) {
    EXPECT_TRUE(true);
}

TEST(FunctionalTest, correct_equilibrium) {
    // check weather or not the pen and computer agree
    handle_t h = 1;
    int dimension = 2;
    int channels = 9;
    array_t pos;
    pos.resize(3);
    pos << 1,2,4;
    auto n = new node(h,dimension,channels,pos,NO_BOUNDARY);
    n->rho = 1;
    n->u.setZero();
    n->population_even= equilibrium(n);
    // should be standard init values
    EXPECT_NEAR(4.0/9,n->population_even(0),1e-10);
    EXPECT_NEAR(1.0/9,n->population_even(1),1e-10);
    EXPECT_NEAR(1.0/9,n->population_even(2),1e-10);
    EXPECT_NEAR(1.0/9,n->population_even(3),1e-10);
    EXPECT_NEAR(1.0/9,n->population_even(4),1e-10);
    EXPECT_NEAR(1.0/36,n->population_even(5),1e-10);
    EXPECT_NEAR(1.0/36,n->population_even(6),1e-10);
    EXPECT_NEAR(1.0/36,n->population_even(7),1e-10);
    EXPECT_NEAR(1.0/36,n->population_even(8),1e-10);
    // other values
    n->rho = 1;
    n->u.setOnes();
    n->population_even= equilibrium(n);
    EXPECT_NEAR(-8.0/9,n->population_even(0),1e-10);
    EXPECT_NEAR(11.0/18,n->population_even(1),1e-10);
    EXPECT_NEAR(11.0/18,n->population_even(2),1e-10);
    EXPECT_NEAR(-1.0/18,n->population_even(3),1e-10);
    EXPECT_NEAR(-1.0/18,n->population_even(4),1e-10);
    EXPECT_NEAR(22.0/36,n->population_even(5),1e-10);
    EXPECT_NEAR(-2.0/36,n->population_even(6),1e-10);
    EXPECT_NEAR(10.0/36,n->population_even(7),1e-10);
    EXPECT_NEAR(-2.0/36,n->population_even(8),1e-10);
    // -1 velocities
    n->rho = 1;
    n->u.setOnes();
    n->u *= -1;
    n->population_even= equilibrium(n);
    EXPECT_NEAR(-8.0/9,n->population_even(0),1e-10);
    EXPECT_NEAR(-1.0/18,n->population_even(1),1e-10);
    EXPECT_NEAR(-1.0/18,n->population_even(2),1e-10);
    EXPECT_NEAR(11.0/18,n->population_even(3),1e-10);
    EXPECT_NEAR(11.0/18,n->population_even(4),1e-10);
    EXPECT_NEAR(10.0/36,n->population_even(5),1e-10);
    EXPECT_NEAR(-2.0/36,n->population_even(6),1e-10);
    EXPECT_NEAR(22.0/36,n->population_even(7),1e-10);
    EXPECT_NEAR(-2.0/36,n->population_even(8),1e-10);
}

TEST(FunctionalTest,correct_macro) {
    // values should just read back
    handle_t h = 1;
    int dimension = 2;
    int channels = 9;
    array_t pos;
    pos.resize(3);
    pos << 1,2,4;
    auto n = new node(h,dimension,channels,pos,NO_BOUNDARY);
    n->rho = 1;
    n->u.setZero();
    n->population_even= equilibrium(n);
    macro(n);
    EXPECT_NEAR(n->rho, 1,1e-10);
    EXPECT_NEAR(n->u(0), 0,1e-10);
    EXPECT_NEAR(n->u(1), 0,1e-10);
    // next
    n->rho = 1;
    n->u.setOnes();
    n->population_even= equilibrium(n);
    macro(n);
    EXPECT_NEAR(n->rho, 1, 1e-10);
    EXPECT_NEAR(n->u(0), 1, 1e-10);
    EXPECT_NEAR(n->u(1), 1, 1e-10);
    // next
    n->rho = 1;
    n->u.setOnes();
    n->u *= 10;
    n->population_even= equilibrium(n);
    macro(n);
    EXPECT_NEAR(n->rho, 1, 1e-10);
    EXPECT_NEAR(n->u(0), 10, 1e-10);
    EXPECT_NEAR(n->u(1), 10, 1e-10);
}

TEST(FunctionalTest, equilibrium123) {
    /// testing the equilbirum function
    handle_t h = 1;
    double rho = 4;
    double ux = 6;
    double uy = 3;
    int dimension = 2;
    int channels = 9;
    array_t pos;
    pos.resize(3);
    pos << 1,2,4;
    auto n = new node(h,dimension,channels,pos,NO_BOUNDARY);
    n->rho = rho;
    n->u << ux,uy;
    n->population_even= equilibrium(n);
    /// expections
    EXPECT_EQ(n->population_even(0), rho * 2/9 *(2 - 3* (ux*ux +uy*uy)));
    EXPECT_EQ(n->population_even(1),rho * 1/18 * (2 + 6*ux + 9*ux*ux - 3*(ux*ux + uy*uy)));
    EXPECT_EQ(n->population_even(2),rho * 1/18 * (2 + 6*uy + 9*uy*uy - 3*(ux*ux + uy*uy)));
    EXPECT_EQ(n->population_even(3),rho * 1/18 * (2 - 6*ux + 9*ux*ux - 3*(ux*ux + uy*uy)));
    EXPECT_EQ(n->population_even(4),rho * 1/18 * (2 - 6*uy + 9*uy*uy - 3*(ux*ux + uy*uy)));
    EXPECT_EQ(n->population_even(5),rho * 1/36 *(1 + 3 *(ux + uy) + 9*ux*uy + 3*(ux*ux + uy*uy)));
    EXPECT_EQ(n->population_even(6),rho * 1/36 *(1 - 3 *(ux - uy) - 9*ux*uy + 3*(ux*ux + uy*uy)));
    EXPECT_EQ(n->population_even(7),rho * 1/36 *(1 - 3 *(ux + uy) + 9*ux*uy + 3*(ux*ux + uy*uy)));
    EXPECT_EQ(n->population_even(8),rho * 1/36 *(1 + 3 *(ux - uy) - 9*ux*uy + 3*(ux*ux + uy*uy)));
}

TEST(FunctionalTest, macro123) {
    handle_t h = 1;
    int dimension = 2;
    int channels = 9;
    array_t pos;
    pos.resize(3);
    pos << 1,2,4;
    auto n = new node(h,dimension,channels,pos,NO_BOUNDARY);
    n->population_even<< 1,2,3,4,5,6,7,8,9;
    macro(n);
    EXPECT_EQ(45, n->rho);
    EXPECT_NEAR(-2.0/45,n->u(0),1e-10);
    EXPECT_NEAR(-6.0/45,n->u(1), 1e-10);
}

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
    simulation sm(&boundaries);
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

TEST(StreamTests, one_D_streaming_channel_three) {
    int size = 20;
    point_t start = {0,0};
    point_t end = {size,0};
    point_t sim_area = {size,1};
    boundaryPointConstructor boundaries(sim_area);
    boundaries.init_structure();
    boundaries.set_point(&start,NO_BOUNDARY);
    boundaries.set_point(&end  ,NO_BOUNDARY);
    simulation sm(&boundaries);
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
    simulation sm(&boundaries,&gen);
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
    simulation sm(&boundaries,&gen);
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
    simulation sm(&boundaries,&gen);
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
    simulation sm(&boundaries,&gen);
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
    simulation sm(&boundaries,&gen);
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
    simulation sm(&boundaries,&gen);
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

TEST(StreamTests, combinded_test_boundary_consistent) {
    // checks leakages form the boundary nodes
    // setup sim
    int size = 5;
    int sim_size = size-2;
    int steps = 10;
    point_t p = {size,size};
    boundaryPointConstructor boundaries(p);
    boundaries.init_quader();
    simulation sim(&boundaries);
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
    /// zero the nodes and then put the codes in
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
    simulation sm(&boundaries,&gen);
    sm.init();
    // zero the data
    for(auto node : sm.nodes) {
        node->population_even.setZero();
        node->population_odd.setZero();
        /** useful debug fragement
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
    simulation sm(&boundaries,&gen);
    sm.init();
    // zero the data
    for(auto node : sm.nodes) {
        node->population_even.setZero();
        node->population_odd.setZero();
        /** useful debug fragement
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
    simulation sm(&boundaries,&gen);
    sm.init();
    // zero the data
    for(auto node : sm.nodes) {
        node->population_even.setZero();
        node->population_odd.setZero();
        //useful debug fragement
        /**
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

TEST(BounceBackTesting, Oblique_six_eight) {
    // we need to generate two nodes and only really consider the bb
    // generate two points we have two bounary nodes and one wet node
    // in this case we are interested in channel 1 to 3 bb
    int size = 3;
    point_t start = {0,size};
    point_t end = {size,0 };
    point_t sim_area = {size,size};
    // need to be able to prime the boundary point construtor
    vector_t node_generation = {1,-1};
    boundaryPointConstructor boundaries(sim_area);
    boundaries.init_structure();
    boundaries.set_point(&start,BOUNCE_BACK);
    boundaries.set_point(&end  ,BOUNCE_BACK);
    nodeGenerator gen(&boundaries);
    gen.set_discovery_vector(node_generation);
    gen.init();
    simulation sm(&boundaries,&gen);
    sm.init();
    // zero the data
    for(auto node : sm.nodes) {
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
    EXPECT_EQ(sm.nodes.at(0)->population_even(8),1);
    EXPECT_EQ(sm.nodes.at(1)->population_even(6),1);
}

TEST(BounceBackTesting, moving) {
    /// test the actual moving code
    /// if this still doesnt work out ill do poisioulle flow...
    // probable cause leakage into channels that have nothing to do
    // in the corners ?!
    // put a bunch of 1 into the top middle node of a 3x3 simspace
    int size = 5;
    int sim_size = size-2;
    double uw = 0.1;
    point_t p = {size,size};
    boundaryPointConstructor boundaries(p);
    boundaries.init_sliding_lid();
    simulation sim(&boundaries);
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

TEST(StreamTests, channel_0_persistent) {
    // test weather or not values in channel 0 persist throu streaming
    int size = 5;
    int steps = 5;
    point_t p = {size,size};
    boundaryPointConstructor boundaries(p);
    boundaries.init_sliding_lid();
    // simulation init
    simulation sim(&boundaries);
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

TEST(FunctionalTest, read_write) {
    int size = 5;
    point_t p = {size,size};
    boundaryPointConstructor boundaries(p);
    boundaries.init_quader();
    /// run the first node generator
    nodeGenerator nodes(&boundaries);
    nodes.set_redo_save(true, true);
    nodes.init();
    /// run a second one but relay on the result of the first one
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

TEST(Orderingtests, bitInterleaving) {
    EXPECT_EQ(bit_interleaving(2,3), 0xE);
    EXPECT_EQ(bit_interleaving_2d(2,3),0xE);
    EXPECT_EQ(bit_interleaving_3d(1,1,1),0x7);
    EXPECT_EQ(bit_interleaving_3d(3,3,3),0x3F);
    EXPECT_EQ(bit_interleaving_3d(7,7,7),0x1FF);
    // maximums
    EXPECT_EQ(bit_interleaving_3d(0x1FFFFF,0x1FFFFF,0x1FFFFF),0x7FFFFFFFFFFFFFFF);
    EXPECT_EQ(bit_interleaving_3d(0x0,0x0,0x0),0x0);
    EXPECT_EQ(bit_interleaving_3d(0xFFFFFF,0xFFFFFF,0xFFFFFF),0x7FFFFFFFFFFFFFFF);
}

TEST(Orderingtests,bit_in_out) {
    uint32_t x = 2;
    uint32_t y = 5;
    uint32_t z = 6;
    // 3d
    uint64_t o = bit_interleaving_3d(x,y,z);
    EXPECT_EQ(x, bit_extraleaving_3d_x(o));
    EXPECT_EQ(y, bit_extraleaving_3d_y(o));
    EXPECT_EQ(z, bit_extraleaving_3d_z(o));
    // 2d
    o = bit_interleaving_2d(x,y);
    EXPECT_EQ(x, bit_extraleaving_2d_x(o));
    EXPECT_EQ(y, bit_extraleaving_2d_y(o));
}

TEST(Orderingtests, sign_reduce) {
    EXPECT_EQ(reduce_32_2(0),0);
    EXPECT_EQ(reduce_32_2(1),1);
    EXPECT_EQ(reduce_32_2(-1),3);
    EXPECT_EQ(reduce_32_2(-0),0);
    //dks
    EXPECT_EQ(reduce_32_2(-4),0x2);
}

TEST(FunctionalTest, fused_collision) {
    // test against original
    double relaxation_time = 0.5;
    point_t pos = {0,0};
    node node_original(1,velocity_set.rows(),velocity_set.cols(),pos,NO_BOUNDARY);
    node node_fused(1,velocity_set.rows(),velocity_set.cols(),pos,NO_BOUNDARY);
    collision(&node_original,relaxation_time);
    fused_collision(&node_fused, relaxation_time);
    // check against each other
    for(int i = 0; i < velocity_set.cols(); ++i) {
        EXPECT_EQ(node_original.population_even(i),node_fused.population_even(i));
    }
}

TEST(FunctionalTest, fused_macro) {
    // test against original
    double relaxation_time = 0.5;
    point_t pos = {0,0};
    node node_original(1,velocity_set.rows(),velocity_set.cols(),pos,NO_BOUNDARY);
    node node_fused(1,velocity_set.rows(),velocity_set.cols(),pos,NO_BOUNDARY);
    // set good values
    node_original.population_even.setOnes();
    node_fused.population_even.setOnes();
    // run functions
    macro(&node_original);
    fused_macro(&node_fused);
    // check
    EXPECT_EQ(node_original.rho,node_fused.rho);
    EXPECT_EQ(node_original.u(0),node_fused.u(0));
    EXPECT_EQ(node_original.u(1),node_fused.u(1));
}

TEST(FunctionalTest, fused_streaming_13) {
    // sizes matters
    int size = 4;
    point_t start = {0,0};
    point_t end = {size,0};
    point_t sim_area = {size,size-1};
    // init
    boundaryPointConstructor boundaries(sim_area);
    boundaries.init_quader();
    boundaries.visualize_2D_boundary();
    nodeGenerator gen(&boundaries);
    gen.init_fused(size);
    gen.visualize_2D_nodes();
    EXPECT_EQ(gen.node_infos.size(),2);
    // gen.visualize_2D_nodes(4);
    simulation sm(&boundaries, &gen);
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

TEST(FunctionalTest, fused_streaming_24) {
    // sizes matters
    int size = 4;
    point_t start = {0,0};
    point_t end = {size,0};
    point_t sim_area = {size-1,size};
    // init
    boundaryPointConstructor boundaries(sim_area);
    boundaries.init_quader();
    //boundaries.visualize_2D_boundary(size);
    nodeGenerator gen(&boundaries);
    gen.init_fused(size);
    //gen.visualize_2D_nodes(4);
    simulation sm(&boundaries, &gen);
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

TEST(FunctionalTest, one_step_macro_collison) {
    // test against original
    double relaxation_time = 0.5;
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
    one_step_macro_collision(&s_node,relaxation_time);
    // compare values
    int o= s_node.offset;
    for(int i = 0; i < velocity_set.cols(); ++i) {
        EXPECT_EQ(node_original.population_even(i),node_fused.population_even(i));
        EXPECT_EQ(s_node.populations(o + i),node_original.population_even(i));
    }
}

TEST(FunctionalTest, oSimu_streaming_13) {
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
    oSimu sm(&boundaries, &gen);
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
    for(auto n : sm.nodes) {
        n->offset = (step & 0x1) * 9;
        sm.streaming(n);
    }
    step++;
    EXPECT_EQ(sm.nodes.at(0)->populations(1 + sm.offset_sim), 1);
    EXPECT_EQ(sm.nodes.at(1)->populations(3 + sm.offset_sim), 1);
    // propagation test
    sm.offset_sim = ((step +1) & 0x1) * 9;
    for(auto n : sm.nodes) {
        n->offset = (step & 0x1) * 9;
        sm.streaming(n);
    }
    step++;
    EXPECT_EQ(sm.nodes.at(0)->populations(3 + sm.offset_sim), 1);
    EXPECT_EQ(sm.nodes.at(1)->populations(1 + sm.offset_sim), 1);
}

TEST(FunctionalTest, oSimu_streaming_24) {
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
    oSimu sm(&boundaries, &gen);
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
    for(auto n : sm.nodes) {
        n->offset = (step & 0x1) * 9;
        sm.streaming(n);
    }
    step++;
    EXPECT_EQ(sm.nodes.at(0)->populations(2 + sm.offset_sim), 1);
    EXPECT_EQ(sm.nodes.at(1)->populations(4 + sm.offset_sim), 1);
    // propagation test
    sm.offset_sim = ((step +1) & 0x1) * 9;
    for(auto n : sm.nodes) {
        n->offset = (step & 0x1) * 9;
        sm.streaming(n);
    }
    step++;
    EXPECT_EQ(sm.nodes.at(0)->populations(4 + sm.offset_sim), 1);
    EXPECT_EQ(sm.nodes.at(1)->populations(2 + sm.offset_sim), 1);
}

TEST(FunctionalTest, oSimu_streaming_57) {
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
    // boundaries.visualize_2D_boundary(size);
    // node init
    nodeGenerator gen(&boundaries);
    gen.init_fused(size);
    // gen.visualize_2D_nodes(size);
    // tests
    oSimu sm(&boundaries, &gen);
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
    for(auto n : sm.nodes) {
        n->offset = (step & 0x1) * 9;
        sm.streaming(n);
    }
    step++;
    EXPECT_EQ(sm.nodes.at(0)->populations(5+ sm.offset_sim), 1);
    EXPECT_EQ(sm.nodes.at(1)->populations(7 + sm.offset_sim), 1);
    // propagation test
    sm.offset_sim = ((step +1) & 0x1) * 9;
    for(auto n : sm.nodes) {
        n->offset = (step & 0x1) * 9;
        sm.streaming(n);
    }
    step++;
    EXPECT_EQ(sm.nodes.at(0)->populations(7 + sm.offset_sim), 1);
    EXPECT_EQ(sm.nodes.at(1)->populations(5 + sm.offset_sim), 1);
}

TEST(FunctionalTest, oSimu_streaming_68) {
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
    boundaries.visualize_2D_boundary();
    // node init
    // first found points seems to be bugged 2.5,2.5
    nodeGenerator gen(&boundaries);
    gen.init_fused(size);
    // gen.visualize_2D_nodes(size);
    oSimu sm(&boundaries, &gen);
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
    for(auto n : sm.nodes) {
        n->offset = (step & 0x1) * 9;
        sm.streaming(n);
    }
    step++;
    EXPECT_EQ(sm.nodes.at(0)->populations(6+ sm.offset_sim), 1);
    EXPECT_EQ(sm.nodes.at(1)->populations(8 + sm.offset_sim), 1);
    // propagation test
    sm.offset_sim = ((step +1) & 0x1) * 9;
    for(auto n : sm.nodes) {
        n->offset = (step & 0x1) * 9;
        sm.streaming(n);
    }
    step++;
    EXPECT_EQ(sm.nodes.at(0)->populations(8 + sm.offset_sim), 1);
    EXPECT_EQ(sm.nodes.at(1)->populations(6 + sm.offset_sim), 1);
}

TEST(FunctionalTest, key_search_functionality) {
    point_t k = {8,8};
    point_t l = {2,3};
    point_t m = {9,2};
    point_t n = {0,0};
    // emplace into the hash map
    pointKeyHash pkh;
    pkh.fill_key(1,k);
    pkh.fill_key(2,l);
    pkh.fill_key(3,m);
    // now find the keys and 0
    EXPECT_EQ(pkh.key_translation(k),1);
    EXPECT_EQ(pkh.key_translation(l),2);
    EXPECT_EQ(pkh.key_translation(m),3);
    EXPECT_EQ(pkh.key_translation(n),0);
}

TEST(FunctionalTest, right_link_number) {
    int step = 0;
    unsigned int size = 4;
    point_t c = {size,size};
    boundaryPointConstructor boundaries(c);
    // boundaries.init_sliding_lid_side_chopped({20,10},30);
    boundaries.init_poiseuille_flow();
    EXPECT_EQ(boundaries.total_boundary_nodes(),(size-1)*4);
    nodeGenerator gen(&boundaries);
    gen.init_fused(size);
    EXPECT_EQ(gen.node_infos.size(), 2*4);
    // check weather or not the links make sense
    int max_handle = 2*4;
    for(auto ni : gen.node_infos) {
        EXPECT_EQ(ni->links.size(),8);
        for(auto l : ni->links) {
            EXPECT_GE(max_handle,l.handle);
            if(max_handle<l.handle) {
                std::cout << ni->handle << std::endl;
                std::cout << ni->position << std::endl;
                std::cout << l.channel << std::endl;
            }
        }
    }
}

TEST(FunctionalTest, periodics_full) {
    // todo implement me
    // todo pressure periodic = periodic right now
    int step = 0;
    unsigned int size = 4;
    point_t c = {size,size};
    boundaryPointConstructor boundaries(c);
    // boundaries.init_sliding_lid_side_chopped({20,10},30);
    boundaries.init_poiseuille_flow();
    EXPECT_EQ(boundaries.total_boundary_nodes(),(size-1)*4);
    nodeGenerator gen(&boundaries);
    gen.init_fused(size);
    EXPECT_EQ(gen.node_infos.size(), 2*4);
    oSimu sim(&boundaries, &gen);
    sim.init();
    // zero all the data and set some in the last two rows and observe where it goes
    for(auto node : sim.nodes) {
        node->populations.setZero();
        // set continous to 4 if periodics work they just ceep being 4
        node->populations(1) = 4;
        node->populations(3) = 3;
    }
    // manual streaming step begin
    sim.offset_sim = ((step +1) & 0x1) * 9;
    for(auto n : sim.nodes) {
        n->offset = (step & 0x1) * 9;
        sim.streaming(n);
    }
    step++;
    // manual streaming step end
    // test
    for (auto n : sim.nodes) {
        EXPECT_EQ(n->populations(0 + sim.offset_sim), 0);
        EXPECT_EQ(n->populations(1 + sim.offset_sim), 4);
        EXPECT_EQ(n->populations(2 + sim.offset_sim), 0);
        EXPECT_EQ(n->populations(3 + sim.offset_sim), 3);
        EXPECT_EQ(n->populations(4 + sim.offset_sim), 0);
        EXPECT_EQ(n->populations(5 + sim.offset_sim), 0);
        EXPECT_EQ(n->populations(6 + sim.offset_sim), 0);
        EXPECT_EQ(n->populations(7 + sim.offset_sim), 0);
        EXPECT_EQ(n->populations(8 + sim.offset_sim), 0);
    }
}

TEST(FunctionalTest, path_magic) {
    auto b = get_base_path();
    EXPECT_GT(sizeof(b),0);
}

TEST(FunctionalTest, bmp_read_32b) {
    // tests general functionality of the bmp image reading compability of the image converter
    auto bmp_32_test_image = get_base_path();
    bmp_32_test_image.append("tests");
    bmp_32_test_image.append("test_images");
    bmp_32_test_image.append("test_32_bit.bmp");
    imageConverter ic(bmp_32_test_image);
    ic.init();
    // check the correct bmp image size should be 800x600x32
    // set 32 bit generics vs 8bit struct
    EXPECT_EQ(ic.bmp.data.size(),800*600*32/8);
    EXPECT_EQ(ic.return_number_of_colors(),205);
}

TEST(FunctionalTest, bmp_read_24b) {
    // test image setup
    auto bmp_24_test_image = get_base_path();
    bmp_24_test_image.append("tests");
    bmp_24_test_image.append("test_images");
    bmp_24_test_image.append("test_24_bit.bmp");
    imageConverter ic(bmp_24_test_image);
    ic.init();
    EXPECT_EQ(ic.bmp.data.size(),800*400*24/8);
    EXPECT_GT(ic.return_number_of_colors(),2);
}

TEST(FunctionalTest, boarder_create_image_raw) {
    // test image setup
    auto bmp_24_test_image = get_base_path();
    bmp_24_test_image.append("tests");
    bmp_24_test_image.append("test_images");
    bmp_24_test_image.append("test_small.bmp");
    imageConverter ic(bmp_24_test_image);
    //
    ic.init();
    ic.raw_run();
    EXPECT_EQ(ic.raw->raw_boundary_points.size(),ic.raw->reformed_boundary_points.size());
    // number should be constant between runs
    EXPECT_EQ(ic.raw->raw_boundary_points.size(), 36);
    EXPECT_EQ(ic.raw->reformed_boundary_points.size(), 36);
    // test correct handle placement for the reformed boundary points
    handle_t start = 0;
    for(auto reformed_bp : ic.raw->reformed_boundary_points) {
        EXPECT_EQ(reformed_bp->h, ++start);
    }
    ic.raw_cleanup();
    // there should be nothing left in the raw data
    EXPECT_EQ(ic.raw->raw_boundary_points.size(), 0);
    EXPECT_EQ(ic.raw->reformed_boundary_points.size(), 0);
    // the b-struct should have size 36
    EXPECT_EQ(ic.boundaries->boundary_structures[0]->boundary_points.size(), 36);
    // from here it should be smooth sailing
    nodeGenerator gen(ic.boundaries);
    // check if the bps are in order
    start = 0;
    int fails = 0;
    for(auto bp : ic.boundaries->boundary_structures[0]->boundary_points) {
        EXPECT_EQ(bp->h, ++start);
        if(bp->h != start) {
            ++fails;
        }
    }
    EXPECT_EQ(fails, 0);
    // ic.boundaries->visualize_2D_boundary(10);
    // fails in check nodes
    gen.init_fused(ic.return_basic_size());
    // gen.visualize_2D_nodes(10);
    EXPECT_EQ(gen.node_infos.size(), 8*8);
}

TEST(FunctionalTest, boarder_create_more_than_one) {
    // test image setup
    auto bmp_24_test_image = get_base_path();
    bmp_24_test_image.append("tests");
    bmp_24_test_image.append("test_images");
    bmp_24_test_image.append("test_image.bmp");
    imageConverter ic(bmp_24_test_image);
    int max_nodes = 80*80; // bad intel joke incuming
    // run the image converter
    ic.init();
    ic.raw_run();
    // check if raw is bigger than reformed
    EXPECT_LT(ic.raw->raw_boundary_points.size(), max_nodes);
    EXPECT_LT(ic.raw->reformed_boundary_points.size(),
              ic.raw->raw_boundary_points.size());
    // run the into boundary_structures reformer
    ic.raw_cleanup();
    // look at the boundary structure visually
    // ic.boundaries->visualize_2D_boundary();
}

TEST(FunctionalTest, raw_reduce_test) {
    // test image setup
    auto bmp_24_test_image = get_base_path();
    bmp_24_test_image.append("tests");
    bmp_24_test_image.append("test_images");
    bmp_24_test_image.append("black.bmp");
    imageConverter ic(bmp_24_test_image);
    // run the initilizers
    ic.init();
    ic.raw_run();
    // reduce should leave nothing for the reformed boundaries
    EXPECT_EQ(ic.raw->raw_boundary_points.size(),80*80);
    EXPECT_EQ(ic.raw->reformed_boundary_points.size(),0);
    // also checks color
    EXPECT_FALSE(ic.check_for_white_wet_nodes());
}

TEST(FunctionalTest, image_u_image_90_right) {
    // boundary to image +90 right now a feature
    // also prob smarter to do with an l to determine filp
    EXPECT_TRUE(true);
    // test image setup
    auto bmp_24_test_image = get_base_path();
    bmp_24_test_image.append("tests");
    bmp_24_test_image.append("test_images");
    bmp_24_test_image.append("test_u.bmp");
    imageConverter ic(bmp_24_test_image);
    // run the functions
    ic.init();
    ic.run();
    // ic.boundaries->visualize_2D_boundary();
}

TEST(FunctionalTest, intersetion_calculation) {
    // std use case
    straight_t ray;
    straight_t surface;
    ray.point = {0,0};
    ray.direction = {1,1};
    surface.point = {1,1};
    surface.direction = {0,1};
    EXPECT_EQ(calculate_intersection(&ray,&surface),1);
    // orthogonal surface to test out
    ray.point = {0,0};
    ray.direction = {1,1};
    surface.point = {1,1};
    surface.direction = {1,-1};
    EXPECT_TRUE(std::isnan(calculate_intersection(&ray,&surface)));
}

TEST(FunctionalTest, image_outer_inner) {
    // test image setupa
    auto bmp_24_test_image = get_base_path();
    bmp_24_test_image.append("tests");
    bmp_24_test_image.append("test_images");
    bmp_24_test_image.append("test_inner_outer.bmp");
    imageConverter ic(bmp_24_test_image);
    // run the functions
    ic.init();
    ic.run();
    ic.boundaries->visualize_2D_boundary();
    nodeGenerator gen(ic.boundaries);
    gen.init_fused(ic.return_basic_size());
    gen.visualize_2D_nodes();
}

// todo look up book boy Wolf Gladrow on forcing term in LB
// todo update fused init to also work with boundary limits

/**
 * algorithm updates for straight.
 * prob best to generate all possible straight lines in a boundary structure
 * calculate intersections with the mass center
 * discard doubles
 * inteseting part delete inner intersections if tere is an outer intersection
 * order?! not sure if that is nessessary
 */

