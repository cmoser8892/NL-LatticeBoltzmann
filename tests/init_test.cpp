#include "simulation.h"
#include "helper_functions.h"
#include "functions.h"
#include "straight.h"
#include <gtest/gtest.h>

TEST(InitTests, basicBoundaryPoints) {
    int size = 5;
    point_t p = {size,size};
    boundaryPointConstructor boundaries(p);
    boundaries.init_quader();
    // should be (size-1)^ 2
    EXPECT_EQ(boundaries.boundary_points.size(), (size-1)*(size-1));
    // check first and last element too
    point_t checker = {0,0};
    EXPECT_EQ(boundaries.boundary_points.at(0)->point,checker);
    checker = {0,1};
    EXPECT_EQ(boundaries.boundary_points.at(boundaries.boundary_points.size()-1)->point,checker);
}

TEST(InitTests, checkNodeGeneration) {
    int size = 5;
    point_t p = {size,size};
    boundaryPointConstructor boundaries(p);
    boundaries.init_quader();
    ///
    nodeGenerator nodes(&boundaries);
    nodes.init();
    // check for right amout of nodes
    EXPECT_EQ(nodes.node_infos.size(),size*size );
}

TEST(InitTests, basicNeighbors) {
    int size = 5;
    point_t p = {size,size};
    boundaryPointConstructor boundaries(p);
    boundaries.init_quader();
    ///
    nodeGenerator nodes(&boundaries);
    nodes.init();
    for(auto n : nodes.node_infos) {
        if(n->type == WET) {
            EXPECT_EQ(n->links.size(),8);
        }
        else if(n->type == DRY) {
            EXPECT_GE(n->links.size(),1);
        }
        else {
            EXPECT_TRUE(false); // impossible condition
        }
    }
}

TEST(InitTests, init_sliding_Lid_boundaries) {
    int size = 5;
    point_t p = {size,size};
    boundaryPointConstructor boundaries(p);
    boundaries.init_sliding_lid();
    //  check weather or not the top nodes are flagged as moving
    for( auto b: boundaries.boundary_points) {
        if(b->point.y() == size -1) {
            EXPECT_EQ(b->type, BOUNCE_BACK_MOVING);
        }
    }
}

TEST(InitTests, init_sliding_Lid_simulation) {
    int size = 5;
    point_t p = {size,size};
    boundaryPointConstructor boundaries(p);
    boundaries.init_sliding_lid();
    simulation sim(&boundaries);
    sim.init();
    //  check weather or not the top nodes are flagged as moving
    for(auto n : sim.nodes) {
        if(n->position(1) == size -1) {
            EXPECT_EQ(n->node_type,DRY);
            EXPECT_EQ(n->boundary_type,BOUNCE_BACK_MOVING);
        }
    }
}

TEST(InitTests, sim_init_correct_values) {
    // when  started velocity should be 0 and rho 1
    int size = 5;
    point_t p = {size,size};
    boundaryPointConstructor boundaries(p);
    boundaries.init_quader();
    //
    simulation sim(&boundaries);
    sim.init();
    // get the data info a flow field
    flowfield_t ux;
    flowfield_t uy;
    flowfield_t rho;
    ux.resize(size,size);
    uy.resize(size,size);
    rho.resize(size,size);
    for(auto node: sim.nodes) {
        // 2 methods that could be made into on, but for some indices
        write_ux(node,&ux);
        write_uy(node,&uy);
        write_rho(node,&rho);
    }
    // check
    for(int i = 0; i < size; ++i) {
        for(int j = 0; j < size; ++j) {
            EXPECT_EQ(0,ux(i,j));
            EXPECT_EQ(0,uy(i,j));
            EXPECT_EQ(1, rho(i,j));
        }
    }
}

TEST(InitTests, boundary_compostion) {
    // checks how many links boundary nodes have corners have 1 (total of 4)
    // adjacent ones 2 (total of 8)
    // the rest has 3
    int size = 12;
    point_t p = {size,size};
    boundaryPointConstructor boundaries(p);
    boundaries.init_sliding_lid();
    /// sims
    simulation sim(&boundaries);
    sim.init();
    // check compostion
    int zeros = 0;
    int ones = 0;
    int twos = 0;
    int threes = 0;
    int totals = 0;
    for (auto node : sim.nodes) {
        if(node->node_type == DRY) {
            totals++;
            int links = node->neighbors.size();
            switch(links) {
            case 1:
                ones++;
                break;
            case 2:
                twos++;
                break;
            case 3:
                threes++;
                break;
            default:
                zeros++;
                break;
            }
        }
        else {
            EXPECT_EQ(node->neighbors.size(),8);
        }
    }
    // checks
    EXPECT_EQ(zeros, 0);
    EXPECT_EQ(ones, 4);
    EXPECT_EQ(twos, 8);
    EXPECT_EQ(threes, 4*(size-1) - 12);
}


TEST(InitTests,simulation_init_run) {
    // init the boundary points -> then init the simulation ( and with it the node generator)
    // this is just a run test to see weather or not anything crashes durin 1 run
    int size = 10;
    int sim_size = size-2;
    int steps = 5;
    point_t p = {size,size};
    boundaryPointConstructor boundaries(p);
    boundaries.init_quader();
    //
    simulation sim(&boundaries);
    sim.init();
    // debug_node(sim.nodes.at(0),true);
    for(int i = 0; i < steps; ++i) {
        // debug one node flow?!
        sim.run();
        // everything should go out and than in...
        for(int i = 0; i < sim_size; ++i) {
            EXPECT_EQ(sim.nodes.at(i)->rho, 1);
        }
    }
}

TEST(InitTests,simulation_sliding_lid_recheck_boundary_flags) {
    int size = 12;
    point_t p = {size,size};
    boundaryPointConstructor boundaries(p);
    boundaries.init_sliding_lid();
    //
    simulation sim(&boundaries);
    sim.init();
    // auxiliary test
    // retest weather or not the top nodes have the right identification
    for( auto n : sim.nodes) {
        if(n->node_type == DRY) {
            if(n->position(1) == size-1) {
                EXPECT_EQ(n->boundary_type, BOUNCE_BACK_MOVING);
            }
            else {
                EXPECT_EQ(n->boundary_type, BOUNCE_BACK);
            }
        }
    }
}

TEST(InitTests, simulation_link_positions) {
    /// check weather or not the links point to the correct node
    int size = 8;
    point_t p = {size,size};
    boundaryPointConstructor boundaries(p);
    boundaries.init_sliding_lid();
    //
    simulation sim(&boundaries);
    sim.init();
    for(auto node : sim.nodes) {
        for(auto link : node->neighbors) {
            point_t position = node->position;
            vector_t vector = velocity_set.col(link.channel);
            point_t expected_position = position + vector;
            unsigned long partner_array_position = link.handle -1;
            point_t partner_position = sim.nodes.at(partner_array_position)->position;
            EXPECT_EQ(partner_position,expected_position);
        }
    }
}

TEST(InitTests, sufaces) {
    unsigned int sub_size = 8;
    point_t p = {sub_size,sub_size};
    boundaryPointConstructor boundaries(p);
    boundaries.init_quader({0,1});
    straightGenerator st(&boundaries);
    st.init();
    // test size
    EXPECT_EQ(boundaries.boundary_points.size(),4*(sub_size-1));
    EXPECT_EQ(st.surfaces.size(),4*(sub_size-1));
    EXPECT_EQ(st.surfaces.size(),boundaries.boundary_points.size());
}

TEST(InitTests, reduced_surface) {
    unsigned int size = 6;
    unsigned int sub_size = 4;
    point_t p = {sub_size,sub_size};
    boundaryPointConstructor boundaries(p);
    boundaries.init_chopped_sliding_lid({1,1},0);
    nodeGenerator gen(&boundaries);
    gen.init(size);
    // chekc if 8x8
    EXPECT_EQ(gen.node_infos.size(), sub_size*sub_size);
    int bb_moving_counter = 0;
    for (auto n : gen.node_infos) {
        if(n->boundary == BOUNCE_BACK_MOVING) {
            bb_moving_counter++;
        }
    }
    EXPECT_EQ(bb_moving_counter,4);
    for(auto n : gen.node_infos) {
        if(n->type == WET) {
            EXPECT_EQ(n->links.size(),8);
        }
    }
}

TEST(InitTests, chopped_boundaries) {
    unsigned int size = 6;
    unsigned int sub_size = 4;
    point_t p = {sub_size,sub_size};
    boundaryPointConstructor boundaries(p);
    boundaries.init_chopped_sliding_lid({1,1},4);
    EXPECT_EQ(boundaries.boundary_points.size(), (sub_size-1)*4);
    // check if even with the bulge the sizes are still the same
    EXPECT_EQ(boundaries.boundary_points.size(),(sub_size-1)*4);
    // generator
    nodeGenerator gen(&boundaries);
    gen.init(size);
    // check correct node amount
    EXPECT_EQ(gen.node_infos.size(), sub_size*sub_size-1);
}

TEST(InitTests, outer_inner_quader) {
    unsigned int size = 10;
    unsigned int outer_size = 7;
    unsigned int inner_size = 3;
    point_t c = {size,size};
    point_t p = {outer_size,outer_size};
    boundaryPointConstructor boundaries(p);
    // boundaries.init_sliding_lid_side_chopped({20,10},30);
    boundaries.init_sliding_lid_inner({1,1},{3,3},{inner_size,inner_size});
    EXPECT_EQ(boundaries.boundary_points.size(), ((outer_size-1) + (inner_size-1))*4);
    //nodeGenerator gen(&boundaries);
    //gen.init(size);
    // EXPECT_EQ(gen.node_infos.size(),outer_size*outer_size - inner_size*inner_size);
}

TEST(NeighbourhoodTests, hash_keys) {
    std::unordered_multimap<handle_t,char> map = {{1,'a'},{1,'b'},{1,'d'},{2,'b'}};
    EXPECT_EQ(map.size(),4);
}

TEST(InitTests, board_creation) {
    // tests the boards creation in terms of the correct size
    unsigned int size = 10;
    nodeGenerator n(nullptr);
    n.board_creation(size);
    EXPECT_EQ(size*size, n.node_infos.size());
}

TEST(InitTests, init_odd_middle) {
    // when initializing with uneven sizes the middle gets lost
    int size = 5;
    point_t p = {size,size};
    boundaryPointConstructor boundaries(p);
    boundaries.init_sliding_lid();
    nodeGenerator nodeGenerator(&boundaries);
    nodeGenerator.init(size);
    EXPECT_EQ(size*size,nodeGenerator.node_infos.size());
    for(auto n : nodeGenerator.node_infos) {
        if(n->type == WET) {
            std::cout << n->position << std::endl;
        }
    }
}

TEST(InitTests, fused_init) {
    int size = 4;
    point_t p = {size,size};
    boundaryPointConstructor boundaries(p);
    boundaries.init_sliding_lid();
    nodeGenerator nodeGenerator(&boundaries);
    nodeGenerator.init_fused(size);
    EXPECT_EQ((size-2)*(size-2),nodeGenerator.node_infos.size());
    // check correct links
    // size first
    for(auto n : nodeGenerator.node_infos) {
        EXPECT_EQ(n->links.size(),8);
    }
    // individual links
    for(auto n : nodeGenerator.node_infos) {
        int own = 0;
        for(auto l : n->links) {
            if(l.handle == n->handle) {
                own++;
            }
        }
        EXPECT_EQ(own, 5);
    }
}

TEST(InitTests, inner_outer_neighbour_test) {
    unsigned int size = 302;
    unsigned int sub_size = 202;
    point_t c = {size,size};
    point_t p = {sub_size,sub_size+20};
    boundaryPointConstructor boundaries(p);
    // boundaries.init_sliding_lid_side_chopped({20,10},30);
    boundaries.init_sliding_lid_inner({10,20},{34,45},{49,52});
    nodeGenerator gen(&boundaries);
    gen.init(size);
    for(auto n : gen.node_infos) {
        if(n->type == WET) {
            if(n->links.size() != 8) {
                std::cout << n->position << std::endl;
            }
            EXPECT_EQ(n->links.size(),8);
        }

    }
}