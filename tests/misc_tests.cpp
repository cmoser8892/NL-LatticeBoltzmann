// Include a library file to make sure proper includes are set
#include "simulation.h"
#include "node.h"
#include "functions.h"
#include "helper_functions.h"
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
    n->data = equilibrium(n);
    // should be standard init values
    EXPECT_NEAR(4.0/9,n->data(0),1e-10);
    EXPECT_NEAR(1.0/9,n->data(1),1e-10);
    EXPECT_NEAR(1.0/9,n->data(2),1e-10);
    EXPECT_NEAR(1.0/9,n->data(3),1e-10);
    EXPECT_NEAR(1.0/9,n->data(4),1e-10);
    EXPECT_NEAR(1.0/36,n->data(5),1e-10);
    EXPECT_NEAR(1.0/36,n->data(6),1e-10);
    EXPECT_NEAR(1.0/36,n->data(7),1e-10);
    EXPECT_NEAR(1.0/36,n->data(8),1e-10);
    // other values
    n->rho = 1;
    n->u.setOnes();
    n->data = equilibrium(n);
    EXPECT_NEAR(-8.0/9,n->data(0),1e-10);
    EXPECT_NEAR(11.0/18,n->data(1),1e-10);
    EXPECT_NEAR(11.0/18,n->data(2),1e-10);
    EXPECT_NEAR(-1.0/18,n->data(3),1e-10);
    EXPECT_NEAR(-1.0/18,n->data(4),1e-10);
    EXPECT_NEAR(22.0/36,n->data(5),1e-10);
    EXPECT_NEAR(-2.0/36,n->data(6),1e-10);
    EXPECT_NEAR(10.0/36,n->data(7),1e-10);
    EXPECT_NEAR(-2.0/36,n->data(8),1e-10);
    // -1 velocities
    n->rho = 1;
    n->u.setOnes();
    n->u *= -1;
    n->data = equilibrium(n);
    EXPECT_NEAR(-8.0/9,n->data(0),1e-10);
    EXPECT_NEAR(-1.0/18,n->data(1),1e-10);
    EXPECT_NEAR(-1.0/18,n->data(2),1e-10);
    EXPECT_NEAR(11.0/18,n->data(3),1e-10);
    EXPECT_NEAR(11.0/18,n->data(4),1e-10);
    EXPECT_NEAR(10.0/36,n->data(5),1e-10);
    EXPECT_NEAR(-2.0/36,n->data(6),1e-10);
    EXPECT_NEAR(22.0/36,n->data(7),1e-10);
    EXPECT_NEAR(-2.0/36,n->data(8),1e-10);
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
    n->data = equilibrium(n);
    macro(n);
    EXPECT_NEAR(n->rho, 1,1e-10);
    EXPECT_NEAR(n->u(0), 0,1e-10);
    EXPECT_NEAR(n->u(1), 0,1e-10);
    // next
    n->rho = 1;
    n->u.setOnes();
    n->data = equilibrium(n);
    macro(n);
    EXPECT_NEAR(n->rho, 1, 1e-10);
    EXPECT_NEAR(n->u(0), 1, 1e-10);
    EXPECT_NEAR(n->u(1), 1, 1e-10);
    // next
    n->rho = 1;
    n->u.setOnes();
    n->u *= 10;
    n->data = equilibrium(n);
    macro(n);
    EXPECT_NEAR(n->rho, 1, 1e-10);
    EXPECT_NEAR(n->u(0), 10, 1e-10);
    EXPECT_NEAR(n->u(1), 10, 1e-10);
}

TEST(FunctionalTest, equilibrium123) {
    /// testing the equilbirum function
    // todo maybe introduce some fuzzing ?!
    handle_t h = 1;
    double rho = 1;
    double ux = 2;
    double uy = 3;
    int dimension = 2;
    int channels = 9;
    array_t pos;
    pos.resize(3);
    pos << 1,2,4;
    auto n = new node(h,dimension,channels,pos,NO_BOUNDARY);
    n->rho = rho;
    n->u << ux,uy;
    n->data = equilibrium(n);
    /// expections
    EXPECT_EQ(n->data(0), rho * 2/9 *(2 - 3* (ux*ux +uy*uy)));
    EXPECT_EQ(n->data(1),rho * 1/18 * (2 + 6*ux + 9*ux*ux - 3*(ux*ux + uy*uy)));
    EXPECT_EQ(n->data(2),rho * 1/18 * (2 + 6*uy + 9*uy*uy - 3*(ux*ux + uy*uy)));
    EXPECT_EQ(n->data(3),rho * 1/18 * (2 - 6*ux + 9*ux*ux - 3*(ux*ux + uy*uy)));
    EXPECT_EQ(n->data(4),rho * 1/18 * (2 - 6*uy + 9*uy*uy - 3*(ux*ux + uy*uy)));
    EXPECT_EQ(n->data(5),rho * 1/36 *(1 + 3 *(ux + uy) + 9*ux*uy + 3*(ux*ux + uy*uy)));
    EXPECT_EQ(n->data(6),rho * 1/36 *(1 - 3 *(ux - uy) - 9*ux*uy + 3*(ux*ux + uy*uy)));
    EXPECT_EQ(n->data(7),rho * 1/36 *(1 - 3 *(ux + uy) + 9*ux*uy + 3*(ux*ux + uy*uy)));
    EXPECT_EQ(n->data(8),rho * 1/36 *(1 + 3 *(ux - uy) - 9*ux*uy + 3*(ux*ux + uy*uy)));
}

TEST(FunctionalTest, macro123) {
    handle_t h = 1;
    int dimension = 2;
    int channels = 9;
    array_t pos;
    pos.resize(3);
    pos << 1,2,4;
    auto n = new node(h,dimension,channels,pos,NO_BOUNDARY);
    n->data << 1,2,3,4,5,6,7,8,9;
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
    boundaries.set_point(&start,NO_BOUNDARY);
    boundaries.set_point(&end  ,NO_BOUNDARY);
    simulation sm(&boundaries);
    sm.init();
    // zero the data
    for(auto node : sm.nodes) {
        node->data.setZero();
        node->copy.setZero();
    }
    //manual streaming channel 1 start at handle 1 at 1,0 to handle 4 at 4,0
    sm.nodes.at(0)->data(1) = 1;
    for(int i = 1; i < size; ++i) {
        // go check
        for(auto node : sm.nodes) {
            int expect = 0;
            if(node->handle == i) {
                expect = 1;
            }
            EXPECT_EQ(expect, node->data(1));
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
    boundaries.set_point(&start,NO_BOUNDARY);
    boundaries.set_point(&end  ,NO_BOUNDARY);
    simulation sm(&boundaries);
    sm.init();
    // zero the data
    for(auto node : sm.nodes) {
        node->data.setZero();
        node->copy.setZero();
    }
    //manual streaming channel 1
    // the -3 is wierd, but if u think about it last tow nodes
    // have handle 20 and 19 in array 19 and 18 so last on is 17 in array
    sm.nodes.at(size-3)->data(3) = 1;
    // aka 20 nodes in total numbered 1 to 20 last two a boundary so we ignore those
    for(int i = size-2; i > 0; --i) {
        // go check
        for(auto node : sm.nodes) {
            // std::cout << node->neighbors.size() << std::endl;
            int expect = 0;
            if(node->handle == i) {
                expect = 1;
            }
            EXPECT_EQ(expect, node->data(3));
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
    boundaries.set_point(&start,NO_BOUNDARY);
    boundaries.set_point(&end  ,NO_BOUNDARY);
    nodeGenerator gen(&boundaries);
    gen.set_discovery_vector(node_generation);
    gen.init();
    simulation sm(&boundaries,&gen);
    sm.init();
    // zero the data
    for(auto node : sm.nodes) {
        node->data.setZero();
        node->copy.setZero();
    }
    //manual streaming channel 1
    // the -3 is wierd, but if u think about it last tow nodes
    // have handle 20 and 19 in array 19 and 18 so last on is 17 in array
    sm.nodes.at(0)->data(channel) = 1;
    // aka 20 nodes in total numbered 1 to 20 last two a boundary so we ignore those
    for(int i = 1; i < size; ++i) {
        // go check
        for(auto node : sm.nodes) {
            // std::cout << node->neighbors.size() << std::endl;
            int expect = 0;
            if(node->handle == i) {
                expect = 1;
            }
            EXPECT_EQ(expect, node->data(channel));
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
    boundaries.set_point(&start,NO_BOUNDARY);
    boundaries.set_point(&end  ,NO_BOUNDARY);
    nodeGenerator gen(&boundaries);
    gen.set_discovery_vector(node_generation);
    gen.init();
    simulation sm(&boundaries,&gen);
    sm.init();
    // zero the data
    for(auto node : sm.nodes) {
        node->data.setZero();
        node->copy.setZero();
    }
    //manual streaming channel 1
    // the -3 is wierd, but if u think about it last tow nodes
    // have handle 20 and 19 in array 19 and 18 so last on is 17 in array
    sm.nodes.at(size-3)->data(channel) = 1;
    // aka 20 nodes in total numbered 1 to 20 last two a boundary so we ignore those
    for(int i = size-2; i > 0; --i) {
        // go check
        for(auto node : sm.nodes) {
            // std::cout << node->neighbors.size() << std::endl;
            int expect = 0;
            if(node->handle == i) {
                expect = 1;
            }
            EXPECT_EQ(expect, node->data(channel));
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
    boundaries.set_point(&start,NO_BOUNDARY);
    boundaries.set_point(&end  ,NO_BOUNDARY);
    nodeGenerator gen(&boundaries);
    gen.set_discovery_vector(node_generation);
    gen.init();
    simulation sm(&boundaries,&gen);
    sm.init();
    // zero the data
    for(auto node : sm.nodes) {
        node->data.setZero();
        node->copy.setZero();
    }
    //manual streaming channel 1
    // the -3 is wierd, but if u think about it last tow nodes
    // have handle 20 and 19 in array 19 and 18 so last on is 17 in array
    sm.nodes.at(0)->data(channel) = 1;
    // aka 20 nodes in total numbered 1 to 20 last two a boundary so we ignore those
    for(int i = 1; i < size; ++i) {
        // go check
        for(auto node : sm.nodes) {
            // std::cout << node->neighbors.size() << std::endl;
            int expect = 0;
            if(node->handle == i) {
                expect = 1;
            }
            EXPECT_EQ(expect, node->data(channel));
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
    boundaries.set_point(&start,NO_BOUNDARY);
    boundaries.set_point(&end  ,NO_BOUNDARY);
    nodeGenerator gen(&boundaries);
    gen.set_discovery_vector(node_generation);
    gen.init();
    simulation sm(&boundaries,&gen);
    sm.init();
    // zero the data
    for(auto node : sm.nodes) {
        node->data.setZero();
        node->copy.setZero();
    }
    //manual streaming channel 1
    // the -3 is wierd, but if u think about it last tow nodes
    // have handle 20 and 19 in array 19 and 18 so last on is 17 in array
    sm.nodes.at(size-3)->data(channel) = 1;
    // aka 20 nodes in total numbered 1 to 20 last two a boundary so we ignore those
    for(int i = size-2; i > 0; --i) {
        // go check
        for(auto node : sm.nodes) {
            // std::cout << node->neighbors.size() << std::endl;
            int expect = 0;
            if(node->handle == i) {
                expect = 1;
            }
            EXPECT_EQ(expect, node->data(channel));
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
    boundaries.set_point(&start,NO_BOUNDARY);
    boundaries.set_point(&end  ,NO_BOUNDARY);
    nodeGenerator gen(&boundaries);
    gen.set_discovery_vector(node_generation);
    gen.init();
    simulation sm(&boundaries,&gen);
    sm.init();
    // zero the data
    for(auto node : sm.nodes) {
        node->data.setZero();
        node->copy.setZero();
    }
    //manual streaming channel 1
    // the -3 is wierd, but if u think about it last tow nodes
    // have handle 20 and 19 in array 19 and 18 so last on is 17 in array
    sm.nodes.at(0)->data(channel) = 1;
    // aka 20 nodes in total numbered 1 to 20 last two a boundary so we ignore those
    for(int i = 1; i < size; ++i) {
        // go check
        for(auto node : sm.nodes) {
            // std::cout << node->neighbors.size() << std::endl;
            int expect = 0;
            if(node->handle == i) {
                expect = 1;
            }
            EXPECT_EQ(expect, node->data(channel));
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
    boundaries.set_point(&start,NO_BOUNDARY);
    boundaries.set_point(&end  ,NO_BOUNDARY);
    nodeGenerator gen(&boundaries);
    gen.set_discovery_vector(node_generation);
    gen.init();
    simulation sm(&boundaries,&gen);
    sm.init();
    // zero the data
    for(auto node : sm.nodes) {
        node->data.setZero();
        node->copy.setZero();
    }
    //manual streaming channel 1
    // the -3 is wierd, but if u think about it last tow nodes
    // have handle 20 and 19 in array 19 and 18 so last on is 17 in array
    sm.nodes.at(size-3)->data(channel) = 1;
    // aka 20 nodes in total numbered 1 to 20 last two a boundary so we ignore those
    for(int i = size-2; i > 0; --i) {
        // go check
        for(auto node : sm.nodes) {
            // std::cout << node->neighbors.size() << std::endl;
            int expect = 0;
            if(node->handle == i) {
                expect = 1;
            }
            EXPECT_EQ(expect, node->data(channel));
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
            node->data = 0;
            node->copy = 0;
            macro(node);
        }
        else {
            node->data = 1;
            node->copy = 1;
            macro(node);
        }
    }
    // 4 corners
    for( int i = 0; i < 4; ++i) {
        for(auto node: sim.nodes) {
            if(node_position_comparison(node,&index_corners.at(i))) {
                // increment
                for(int j = 0; j < node->data.size(); ++j) {
                    node->data(j) = index_codes.at(i) * 10 + j;
                }
                node->copy = node->data;
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

TEST(StreamTests, DISABLED_combinded_test_1_step) {
    /// hard to get that one right
    // set some codes in the corners and see how it moves around
    // also checks leakages form the boundary nodes
    // setup sim
    int size = 6;
    int sim_size = size-2;
    int steps = 1;
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
    array_t  corner_2(2); corner_2 << sim_size,1; ;
    int code_107 = 107;
    array_t  corner_3(2); corner_3 << sim_size,sim_size;
    int code_126 = 126;
    array_t  corner_4(2); corner_4 << 1,sim_size;
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
            node->data = 0;
            node->copy = 0;
            macro(node);
        }
        else {
            node->data = 1;
            node->copy = 1;
            macro(node);
        }
    }
    // 4 corners
    for( int i = 0; i < 4; ++i) {
        for(auto node: sim.nodes) {
            if(node_position_comparison(node,&index_corners.at(i))) {
                // increment
                for(int j = 0; j < node->data.size(); ++j) {
                    node->data(j) = index_codes.at(i) * 10 + j;
                }
                node->copy = node->data;
            }
        }
    }
    // run checks rest is preamble
    for(int i = 0; i < steps; i++) {
        sim.streaming_step_1();
        sim.streaming_step_2();
        sim.bounce_back();
        // check
        for(auto node : sim.nodes) {
            debug_node(node,true);
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
    boundaries.set_point(&start,BOUNCE_BACK);
    boundaries.set_point(&end  ,BOUNCE_BACK);
    simulation sm(&boundaries);
    sm.init();
    // zero the data
    for(auto node : sm.nodes) {
        node->data.setZero();
        node->copy.setZero();
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
    sm.nodes.at(0)->data(3) = 1;
    sm.nodes.at(1)->data(1) = 1;
    sm.streaming_step_1();
    sm.streaming_step_2();
    sm.bounce_back();
    EXPECT_EQ(sm.nodes.at(1)->data(3),1);
    EXPECT_EQ(sm.nodes.at(0)->data(1),1);
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
    boundaries.set_point(&start,BOUNCE_BACK);
    boundaries.set_point(&end  ,BOUNCE_BACK);
    nodeGenerator gen(&boundaries);
    gen.set_discovery_vector(node_generation);
    gen.init();
    simulation sm(&boundaries,&gen);
    sm.init();
    // zero the data
    for(auto node : sm.nodes) {
        node->data.setZero();
        node->copy.setZero();
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
    sm.nodes.at(0)->data(2) = 1;
    sm.nodes.at(1)->data(4) = 1;
    sm.streaming_step_1();
    sm.streaming_step_2();
    sm.bounce_back();
    EXPECT_EQ(sm.nodes.at(1)->data(2),1);
    EXPECT_EQ(sm.nodes.at(0)->data(4),1);
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
    boundaries.set_point(&start,BOUNCE_BACK);
    boundaries.set_point(&end  ,BOUNCE_BACK);
    nodeGenerator gen(&boundaries);
    gen.set_discovery_vector(node_generation);
    gen.init();
    simulation sm(&boundaries,&gen);
    sm.init();
    // zero the data
    for(auto node : sm.nodes) {
        node->data.setZero();
        node->copy.setZero();
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
    sm.nodes.at(0)->data(7) = 1;
    sm.nodes.at(1)->data(5) = 1;
    sm.streaming_step_1();
    sm.streaming_step_2();
    sm.bounce_back();
    EXPECT_EQ(sm.nodes.at(0)->data(5),1);
    EXPECT_EQ(sm.nodes.at(1)->data(7),1);
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
    boundaries.set_point(&start,BOUNCE_BACK);
    boundaries.set_point(&end  ,BOUNCE_BACK);
    nodeGenerator gen(&boundaries);
    gen.set_discovery_vector(node_generation);
    gen.init();
    simulation sm(&boundaries,&gen);
    sm.init();
    // zero the data
    for(auto node : sm.nodes) {
        node->data.setZero();
        node->copy.setZero();
    }
    // set the value of channel 1 in the middle node to 1
    // should reappear in channel 3 and vice versa
    // we use two points with handle 1 & 2 to check the behaviour of bb in channel 1 and 3
    sm.nodes.at(0)->data(6) = 1;
    sm.nodes.at(1)->data(8) = 1;
    sm.streaming_step_1();
    sm.streaming_step_2();
    sm.bounce_back();
    EXPECT_EQ(sm.nodes.at(0)->data(8),1);
    EXPECT_EQ(sm.nodes.at(1)->data(6),1);
}

TEST(BounceBackTesting, moving) {
    /// test the actual moving code
    /// if this still doesnt work out ill do poisioulle flow...
    // probable cause leakage into channels that have nothing to do
    // in the corners ?!
    // put a bunch of 1 into the top middle node of a 3x3 simspace
    int size = 5;
    point_t p = {size,size};
    boundaryPointConstructor boundaries(p);
    boundaries.init_sliding_lid();
    simulation sim(&boundaries);
    sim.init();
    // zero the data
    for(auto node : sim.nodes) {
        node->data.setZero();
        node->copy.setZero();
    }
    // points of interest and codes
    std::vector<handle_t> handles;
    std::vector<array_t> index_corners;
    std::vector<int> index_codes;
    array_t  poi_1(2); poi_1 << 2,1;
    int code_104 = 104;
    array_t  poi_2(2); poi_2 << 2,3; ;
    int code_107 = 107;
    index_corners.push_back(poi_1);
    index_corners.push_back(poi_2);
    index_codes.push_back(code_104);
    index_codes.push_back(code_107);
    // place the codes into the nodes
    for( int i = 0; i < 2; ++i) {
        for(auto node: sim.nodes) {
            if(node_position_comparison(node,&index_corners.at(i))) {
                // safe the handle for checking results
                handles.push_back(node->handle);
                // put data in
                for(int j = 0; j < node->data.size(); ++j) {
                    // node->data(j) = index_codes.at(i) * 10 + j;
                }
                node->copy = node->data;
            }
        }
    }
    // run the bb sequence
    sim.streaming_step_1();
    sim.streaming_step_2();
    sim.bounce_back();
    // check results
    int position_poi_1 = int(handles.at(0)-1);
    int position_poi_2 = int(handles.at(1)-1);
    for(auto node : sim.nodes) {
        debug_node(node,true);
    }
}