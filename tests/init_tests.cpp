#include "simulation.h"
#include "two_step_simulation.h"
#include "one_step_simulation.h"
#include "helper_functions.h"
#include "functions.h"
#include "straight.h"
#include <gtest/gtest.h>

/**
 * Tests the init quader method and checks weather or not it has the right amount of nodes.
 * @test
 * @see boundaryPointConstructor::init_quader()
 */
TEST(InitTests, basicBoundaryPoints) {
    int size = 5;
    point_t p = {size,size};
    boundaryPointConstructor boundaries(p);
    boundaries.init_quader();
    // should be (size-1)^ 2
    EXPECT_EQ(boundaries.total_boundary_nodes(), (size-1)*(size-1));
}

/**
 * Checks if the nodeGenerator intis the right amount of nodes.
 * @test
 * @see boundaryPointConstructor::init()
 * @see nodeGnerator::init()
 */
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

/**
 * Checks the amount of links a node has in relation to it's position.
 * @test
 * @see boundaryPointConstructor::init()
 * @see nodeGnerator::init()
 */
TEST(InitTests, basicNeighbors) {
    int size = 5;
    point_t p = {size,size};
    boundaryPointConstructor boundaries(p);
    boundaries.init_quader();
    //
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

/**
 * Checks weather or not sliding lid has the right amount of BOUNCE_BACK_MOVING nodes.
 * @test
 * @see boundaryPointConstructor::init_slinding_lid()
 */
TEST(InitTests, init_sliding_Lid_boundaries) {
    int size = 5;
    point_t p = {size,size};
    boundaryPointConstructor boundaries(p);
    boundaries.init_sliding_lid();
    //  check weather or not the top nodes are flagged as moving
    for(auto bs : boundaries.boundary_structures) {
        for( auto b: bs->boundary_points) {
            if(b->point.y() == size -1) {
                EXPECT_EQ(b->type, BOUNCE_BACK_MOVING);
            }
        }
    }
}

/**
 * Tests weather a sliding lid Simulation is set up with the correct amount of BOUNCE_BACK_MOVING nodes.
 * @test
 * @see basicSimulation::init()
 */
TEST(InitTests, init_sliding_Lid_simulation) {
    int size = 5;
    point_t p = {size,size};
    boundaryPointConstructor boundaries(p);
    boundaries.init_sliding_lid();
    basicSimulation sim(&boundaries);
    sim.init();
    // check weather or not the top nodes are flagged as moving
    for(auto n : sim.nodes) {
        if(n->position(1) == size -1) {
            EXPECT_EQ(n->node_type,DRY);
            EXPECT_EQ(n->boundary_type,BOUNCE_BACK_MOVING);
        }
    }
}

/**
 * Checks weather the intial values for rho, ux and uy are set correctly.
 * @test
 * @see basicSimulation::init()
 */
TEST(InitTests, sim_init_correct_values) {
    // when  started velocity should be 0 and rho 1
    int size = 5;
    point_t p = {size,size};
    boundaryPointConstructor boundaries(p);
    boundaries.init_quader();
    //
    basicSimulation sim(&boundaries);
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

/**
 * Tests the links each boundary node has and tests against expected total amount with that meny links.
 * @test
 * @see basicSimulation::init()
 */
TEST(InitTests, boundary_compostion) {
    // checks how many links boundary nodes have corners have 1 (total of 4)
    // adjacent ones 2 (total of 8)
    // the rest has 3
    int size = 12;
    point_t p = {size,size};
    boundaryPointConstructor boundaries(p);
    boundaries.init_sliding_lid();
    // sims
    basicSimulation sim(&boundaries);
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

/**
 * Tests weather or not anything crashes in the init to run process.
 * We just init a quader so nothing actually happens.
 * @test
 */
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
    basicSimulation sim(&boundaries);
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

/**
 * Tests if all the boundary flags are set correctly from the boundaryPointConstructor to simulation.
 * @test
 */
TEST(InitTests,simulation_sliding_lid_recheck_boundary_flags) {
    int size = 12;
    point_t p = {size,size};
    boundaryPointConstructor boundaries(p);
    boundaries.init_sliding_lid();
    //
    basicSimulation sim(&boundaries);
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

/**
 * Tests weather or not the neighborhood of each node is established correctly.
 * @test
 * @see neighbourhood::determine_neighbors()
 */
TEST(InitTests, simulation_link_positions) {
    // check weather or not the links point to the correct node
    int size = 8;
    point_t p = {size,size};
    boundaryPointConstructor boundaries(p);
    boundaries.init_sliding_lid();
    //
    basicSimulation sim(&boundaries);
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

/**
 * Checks weather the amount of surface created from boundary points is correct.
 * @test
 * @see straightGenerator::init()
 */
TEST(InitTests, sufaces) {
    unsigned int sub_size = 8;
    point_t p = {sub_size,sub_size};
    boundaryPointConstructor boundaries(p);
    boundaries.init_quader();
    // boundaries.visualize_2D_boundary();
    straightGenerator st(&boundaries);
    st.init();
    // test size
    EXPECT_EQ(boundaries.total_boundary_nodes(),4*(sub_size-1));
    EXPECT_EQ(st.surfaces.size(),4);
    // test the lengths of the vectors
    for(auto s : st.surfaces) {
        double length = std::abs(s->max_t)*s->direction.norm();
        EXPECT_EQ(length,7);
    }

}

/**
 * Test the first method that has not straight boundary.
 * @test
 * @see boundaryPointConstructor::init_chopped_sliding_lid
 */
TEST(InitTests, reduced_size) {
    unsigned int size = 6;
    unsigned int sub_size = 4;
    point_t c = {size, size};
    point_t p = {sub_size,sub_size};
    boundaryPointConstructor boundaries(c);
    boundaries.init_chopped_sliding_lid({1,1},p,0);
    // boundaries.visualize_2D_boundary();
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

/**
 * Tests the chopped boundary function.
 * @test
 * @see boundaryPointConstructor::init_chopped_sliding_lid()
 */
TEST(InitTests, chopped_boundaries) {
    unsigned int size = 6;
    unsigned int sub_size = 4;
    point_t c = {size,size};
    point_t p = {sub_size,sub_size};
    boundaryPointConstructor boundaries(c);
    boundaries.init_chopped_sliding_lid({1,1},p,4);
    // boundaries.visualize_2D_boundary();
    EXPECT_EQ(boundaries.total_boundary_nodes(), (sub_size-1)*4);
    // check if even with the bulge the sizes are still the same
    EXPECT_EQ(boundaries.total_boundary_nodes(),(sub_size-1)*4);
    // generator
    nodeGenerator gen(&boundaries);
    gen.init(size);
    // check correct node amount
    EXPECT_EQ(gen.node_infos.size(), sub_size*sub_size-1);
}

/**
 * Tests the amount of boundary points for a sliding lid with a hole.
 * @test
 * @see boundaryPointConstructor::init_sliding_lid_inner()
 */
TEST(InitTests, outer_inner_quader) {
    unsigned int size = 10;
    unsigned int outer_size = 7;
    unsigned int inner_size = 3;
    point_t c = {size,size};
    point_t p = {outer_size,outer_size};
    boundaryPointConstructor boundaries(p);
    // boundaries.init_sliding_lid_side_chopped({20,10},30);
    boundaries.init_sliding_lid_inner({1,1},p,{3,3},{inner_size,inner_size});
    EXPECT_EQ(boundaries.total_boundary_nodes(), ((outer_size-1) + (inner_size-1))*4);
}

/**
 * Checks out how and unordered multimap works.
 * @test
 */
TEST(NeighbourhoodTests, hash_keys) {
    std::unordered_multimap<handle_t,char> map = {{1,'a'},{1,'b'},{1,'d'},{2,'b'}};
    EXPECT_EQ(map.size(),4);
}

/**
 * Tests the board creation from the nodeGenerator
 * @test
 * @see nodeGenerator::board_creation()
 */
TEST(InitTests, board_creation) {
    // tests the boards creation in terms of the correct size
    unsigned int size = 10;
    boundaryPointConstructor* bs = nullptr;
    nodeGenerator n(bs);
    n.board_creation(size);
    EXPECT_EQ(size*size, n.node_infos.size());
}

/**
 * Tests against where one node that is also the mass_center of the straghtGenerator gets lost.
 * @test
 * @see straightGenerator::calculate_number_intersections()
 */
TEST(InitTests, init_odd_middle) {
    // when initializing with uneven sizes the middle gets lost
    int size = 5;
    point_t p = {size,size};
    boundaryPointConstructor boundaries(p);
    boundaries.init_sliding_lid();
    nodeGenerator nodeGenerator(&boundaries);
    nodeGenerator.init(size);
    EXPECT_EQ(size*size,nodeGenerator.node_infos.size());
}

/**
 * Tests the fused init method of the nodeGenerator.
 * @test
 * @note the boundary points get optimized away
 * @see nodeGenerator::init_fused()
 */
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

/**
 * There was a bug where some nodes unexpectedly went missing for a inner outer boundary clouds.
 * Root cause was the connection of boundary points that do not belong to the same structure.
 * Fixed with the introduction of boundaryStructures.
 * @test
 * @see boundaryStructure
 */
TEST(InitTests, inner_outer_neighbour_test) {
    // still kinda spooked will try out the write out of the rho field
    unsigned int size = 30;
    unsigned int sub_size = 20;
    point_t c = {size,size};
    point_t p = {sub_size +5,sub_size};
    point_t k = {4,10};
    boundaryPointConstructor boundaries(c);
    // boundaries.init_sliding_lid_side_chopped({20,10},30);
    boundaries.init_sliding_lid_inner({3,5},p,{5,7},k);
    // boundaries.visualize_2D_boundary();
    nodeGenerator gen(&boundaries);
    gen.init(size);
    // gen.visualize_2D_nodes();
    int expected_total_node_number = p.x()*p.y() - (k.x()-2)*(k.y()-2);
    EXPECT_EQ(expected_total_node_number,gen.node_infos.size());
    int number_nodes = 0;
    int number_wet_nodes = 0;
    int number_dry_nodes = 0;
    int broken_nodes = 0;
    int stalkers = 0;
    for(auto n : gen.node_infos) {
        number_nodes++;
        if(n->type == WET) {
            number_wet_nodes++;
            if(n->links.size() != 8) {
                std::cout << n->position.x() << " " << n->position.y() << "; " << n->links.size() << std::endl << std::endl;
                broken_nodes++;
            }
            //EXPECT_EQ(n->links.size(),8);
        }
        else if(n->type == DRY) {
            number_dry_nodes++;
        }
        else {
            stalkers++;
        }
    }

    // there are 4 nodes on the loose
    EXPECT_EQ(broken_nodes,0);
    EXPECT_EQ(stalkers,0);
    EXPECT_EQ(expected_total_node_number,number_nodes);
    EXPECT_EQ(number_dry_nodes,boundaries.total_boundary_nodes());
    EXPECT_EQ(number_wet_nodes, expected_total_node_number - boundaries.total_boundary_nodes());
    // sanity check
    EXPECT_EQ(number_nodes, number_dry_nodes + number_wet_nodes);
}

/**
 * Tests against the unclose boundaryStructure bug.
 * @test
 * @see boundaryStructure
 */
TEST(InitTests, init_out_inner_rho_writeout) {
    // checks number of rho zeros
    unsigned int size = 30;
    unsigned int sub_size = 20;
    point_t c = {size,size};
    point_t p = {sub_size,sub_size+2};
    point_t k = {5,6};
    boundaryPointConstructor boundaries(p);
    // boundaries.init_sliding_lid_side_chopped({20,10},30);
    boundaries.init_sliding_lid_inner({3,5},p,{9,7},k);
    nodeGenerator gen(&boundaries);
    gen.init(size);
    basicSimulation sim(&boundaries,&gen);
    sim.init();
    flowfield_t rho;
    rho.resize(size,size);
    rho.setZero();
    for(auto n : sim.nodes) {
        write_rho(n,&rho);
    }
    // go through the flow-field
    int expected_total_node_number = p.x()*p.y() - (k.x()-2)*(k.y()-2);
    int zeros = 0;
    int ones = 0;
    int errors = 0;
    for(int i = 0; i < size;++i) {
        for(int j = 0; j < size; ++j) {
            if(rho(i,j) == 0) {
                zeros++;
            }
            else if(rho(i,j) == 1) {
                ones++;
            }
            else {
                errors++;
            }
        }
    }
    EXPECT_EQ(expected_total_node_number,ones);
    EXPECT_EQ(errors,0);
}

/**
 * Another test against the bug with missing boundary structures.
 * @test
 * @note Root causues was the introduction of cross connected boundary structures resulting for just using boundaryPoints.
 * @see boundaryStructure
 */
TEST(InitTests, inner_outer_master_test) {
    unsigned int size = 10;
    unsigned int sub_size = 8;
    unsigned int white_list = 6;
    unsigned int inner_size = 4;
    point_t c = {size,size}; // size or the canvas
    point_t p = {sub_size,sub_size}; // size of the  quader
    point_t k = {inner_size,inner_size}; // size of the inner one
    boundaryPointConstructor boundaries(c);
    // boundaries.init_sliding_lid_side_chopped({20,10},30);
    boundaries.init_sliding_lid_inner({1,1},p,{3,3},k);
    // boundaries.visualize_2D_boundary();
    EXPECT_EQ(boundaries.total_boundary_nodes(),(sub_size-1)*4 + (inner_size -1)*4);
    nodeGenerator gen(&boundaries);
    // check boundaries right at least
    gen.init(size);
    // gen.visualize_2D_nodes();
    // boundary point sanity check
    int dry_nodes_number = 0;
    for(auto node : gen.node_infos) {
        if(node->type == DRY) {
            dry_nodes_number++;
        }
    }
    EXPECT_EQ(boundaries.total_boundary_nodes(),dry_nodes_number);
    // wet node check
    int wet_nodes_number = 0;
    int good_wet_nodes = 0;
    int bad_wet_nodes = 0;
    // miss using boundary nodes to whitelist the correct nodes in the mass to narrow down the error mass
    point_t w = {white_list,white_list};
    boundaryPointConstructor white_list_nodes(c);
    white_list_nodes.init_quader({2,2},w);
    // white_list_nodes.visualize_2D_boundary(size);
    EXPECT_EQ((white_list-1)*4, white_list_nodes.total_boundary_nodes());
    for(auto node : gen.node_infos) {
        if(node->type == WET) {
            wet_nodes_number++;
            // check the whitelist for correct nodes
            bool good_node = false;
            for(auto bs: white_list_nodes.boundary_structures) {
                for(auto white : bs->boundary_points) {
                    if(point_t(node->position) == white->point) {
                        good_wet_nodes++;
                        good_node  = true;
                        break;
                    }
                }
            }
            if(!good_node) {
                bad_wet_nodes++;
                std::cout << node->position.x() << " " << node->position.y() << std::endl << std::endl;
            }
        }
    }
    EXPECT_EQ((white_list-1)*4, good_wet_nodes);
    EXPECT_EQ(0,bad_wet_nodes);
}

/**
 * Checks weather or not the correct boundary labels got set while getting rid of all the bounce back nodes.
 * @test
 * @see nodeGenerator::init_fused()
 */
TEST(InitTests, fused_init_correct_boundary_lables) {
    unsigned int size = 6;
    point_t c = {size,size};
    boundaryPointConstructor boundaries(c);
    boundaries.init_sliding_lid();
    nodeGenerator gen(&boundaries);
    gen.init_fused(size);
    // check correct
    EXPECT_EQ(gen.node_infos.size(), (size-2)*(size-2));
    int no_boundary = 0;
    int bounce_back = 0;
    int bounce_back_moving = 0;
    for(auto n : gen.node_infos) {
        if(n->boundary == NO_BOUNDARY) {
            no_boundary++;
        }
        else if(n->boundary == BOUNCE_BACK) {
            bounce_back++;
        }
        else if(n->boundary == BOUNCE_BACK_MOVING) {
            bounce_back_moving++;
        }
        else {
            // error
            EXPECT_TRUE(false);
        }
    }
    EXPECT_EQ(no_boundary,4);
    EXPECT_EQ(bounce_back,8);
    EXPECT_EQ(bounce_back_moving,4);
}

/**
 * Tests for the correct amount of nodes while using fused init.
 * @test
 * @see nodeGenerator::init_fused()
 */
TEST(InitTests, fused_inner_outer_init_simple) {
    unsigned int size = 10;
    unsigned int sub_size = 8;
    unsigned int white_list = 6;
    unsigned int inner_size = 4;
    point_t c = {size,size}; // size or the canvas
    point_t p = {sub_size,sub_size}; // size of the  quader
    point_t k = {inner_size,inner_size}; // size of the inner one
    boundaryPointConstructor boundaries(c);
    // boundaries.init_sliding_lid_side_chopped({20,10},30);
    boundaries.init_sliding_lid_inner({1,1},p,{3,3},k);
    // boundaries.visualize_2D_boundary(size);
    EXPECT_EQ(boundaries.total_boundary_nodes(),(sub_size-1)*4 + (inner_size -1)*4);
    nodeGenerator gen(&boundaries);
    // check boundaries right at least
    gen.init_fused(size);
    // gen.visualize_2D_nodes(size);
    int expected_total_node_number = (p.x()-2)*(p.y()-2) - (k.x()*k.y());
    EXPECT_EQ(gen.node_infos.size(),expected_total_node_number);
}

/**
 * Tests fused init with an outer and inner structure set onto the canvas.
 * @test
 * @see basicSimulation::fused_init()
 */
TEST(InitTests, fused_inner_outer_init_shifted) {
    // checks number of rho zeros
    unsigned int size = 30;
    unsigned int sub_size = 20;
    point_t c = {size,size};
    point_t p = {sub_size,sub_size+2};
    point_t k = {5,6};
    boundaryPointConstructor boundaries(p);
    // boundaries.init_sliding_lid_side_chopped({20,10},30);
    boundaries.init_sliding_lid_inner({3,5},p,{9,7},k);
    int counter = 0;
    for(auto bs : boundaries.boundary_structures) {
        for(auto bp : bs->boundary_points) {
            if(bp->type == BOUNCE_BACK_MOVING) {
                ++counter;
            }
        }
    }
    EXPECT_EQ(counter, sub_size);
    nodeGenerator gen(&boundaries);
    gen.init_fused(size);
    basicSimulation sim(&boundaries,&gen);
    sim.fused_init();
    // check
    int expected_total_node_number = (p.x()-2)*(p.y()-2) - (k.x()*k.y());
    EXPECT_EQ(sim.nodes.size(),expected_total_node_number);
    // nessessary to check the boundary nodes too
}

/**
 * Checks if the poiseulle flow was set up correctly.
 * @test
 * @see boundaryPointConstructor::init_poiseuille_flow(
 */
TEST(InitTests, p_flow) {
    unsigned int size = 30;
    point_t s = {size,size-10};
    boundaryPointConstructor boundaries(s);
    boundaries.init_poiseuille_flow();
    EXPECT_EQ(boundaries.total_boundary_nodes(),30+30+18+18);
    // check if we got 20+20 p flow nodes
    int pb = 0;
    for(auto node : boundaries.boundary_structures[0]->boundary_points) {
        if(node->type== PRESSURE_PERIODIC) {
            ++pb;
        }
    }
    int wet = 0;
    int dry = 0;
    for(auto node : boundaries.boundary_structures[0]->boundary_points) {
        if(node->dw == DRY) {
            dry++;
        }
        else if(node->dw == WET) {
            wet++;
        }
        else {
            EXPECT_TRUE(false);
        }
    }
    EXPECT_EQ(pb, 18+18);
    EXPECT_EQ(wet, 18+18);
    EXPECT_EQ(dry, 30+30);
}

/**
 * Tests the individual corner creation function in the boundaryPointConstructor.
 * @test
 * @see boundaryPointConstructor::corner_creation()
 */
TEST(InitTests, corners) {
    // test all corners if works
    int size = 3;
    point_t s = {size,size};
    boundaryPointConstructor b(s);
    vector_t dir = {1,1};
    point_t start = {1,1};
    b.init_structure();
    b.corner_creation(dir,&start,BOUNCE_BACK);
    // b.visualize_2D_boundary(3);
    EXPECT_EQ(b.total_boundary_nodes(),2);
    b.delete_structures();
    EXPECT_EQ(b.total_boundary_nodes(),0);
    dir = {-1,1};
    start = {1,1};
    b.init_structure();
    b.corner_creation(dir,&start,BOUNCE_BACK);
    // b.visualize_2D_boundary(3);
    EXPECT_EQ(b.total_boundary_nodes(),2);
    b.delete_structures();
    EXPECT_EQ(b.total_boundary_nodes(),0);
    dir = {-1,-1};
    start = {1,1};
    b.init_structure();
    b.corner_creation(dir,&start,BOUNCE_BACK);
    // b.visualize_2D_boundary(3);
    EXPECT_EQ(b.total_boundary_nodes(),2);
    b.delete_structures();
    EXPECT_EQ(b.total_boundary_nodes(),0);
    dir = {1,-1};
    start = {1,1};
    b.init_structure();
    b.corner_creation(dir,&start,BOUNCE_BACK);
    // b.visualize_2D_boundary(3);
    EXPECT_EQ(b.total_boundary_nodes(),2);
    b.delete_structures();
    EXPECT_EQ(b.total_boundary_nodes(),0);
}

/**
 * Tests weather or not the straight generator can handle unordered boundary points.
 * @test
 * @see straightGenerator::init()
 */
TEST(InitTests, straight_unordered) {
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
    point = {0,1};
    boundaries.set_point(&point,BOUNCE_BACK);
    point = {2,0};
    boundaries.set_point(&point,BOUNCE_BACK);
    point = {0,2};
    boundaries.set_point(&point,BOUNCE_BACK);
    point = {3,0};
    boundaries.set_point(&point,BOUNCE_BACK);
    point = {3,1};
    boundaries.set_point(&point,BOUNCE_BACK);
    point = {0,3};
    boundaries.set_point(&point,BOUNCE_BACK);
    point = {3,2};
    boundaries.set_point(&point,BOUNCE_BACK);
    point = {1,3};
    boundaries.set_point(&point,BOUNCE_BACK);
    point = {2,2};
    boundaries.set_point(&point,BOUNCE_BACK);
    point = {2,3};
    boundaries.set_point(&point,BOUNCE_BACK);
    point = {1,1};
    boundaries.set_point(&point,BOUNCE_BACK);
    // run the straight gen
    straightGenerator straight(&boundaries);
    straight.init();
    // checks number and length of each generated straight
    EXPECT_EQ(boundaries.total_boundary_nodes(),12);
    for(auto lines : straight.surfaces) {
        EXPECT_EQ(lines->direction.norm(),1);
    }
    handle_t expected_handle = 1;
    for(auto bp : boundaries.boundary_structures[0]->boundary_points) {
        EXPECT_EQ(expected_handle, bp->h);
        expected_handle++;
    }
}

/**
 * When we want to handle dry and wet boundaries they have to be ordered so that the wet nodes come first and dont get lost.
 * @test
 * @see boundaryPointConstructor::init_poiseuille_flow()
 */
TEST(InitTests, ordering_boundaries) {
    unsigned int size = 4;
    point_t c = {size,size};
    boundaryPointConstructor boundaries(c);
    // boundaries.init_sliding_lid_side_chopped({20,10},30);
    boundaries.init_poiseuille_flow();
    bool wet_dry_change_around = false;
    // this actually copies the content in the vector over to bs ?!
    auto bs = boundaries.boundary_structures[0]->boundary_points;
    for(auto bp : bs) {
        if(bp->dw == DRY) {
            wet_dry_change_around = true;
        }
        if(wet_dry_change_around) {
            EXPECT_EQ(bp->dw,DRY);
        }
        else {
            EXPECT_EQ(bp->dw, WET);
        }
    }
}

/**
 * Tests weather or not a simulation can be init with wet boundary nodes.
 * Checks if the sizes are ok.
 * @test
 */
TEST(InitTests, special_case_wet_boundary) {
    // when labeling some boundary nodes as wet we get an error in reduce boundary neighborhood
    unsigned int size = 4;
    point_t c = {size,size};
    boundaryPointConstructor boundaries(c);
    // boundaries.init_sliding_lid_side_chopped({20,10},30);
    boundaries.init_poiseuille_flow();
    EXPECT_EQ(boundaries.total_boundary_nodes(),(size-1)*4);
    nodeGenerator gen(&boundaries);
    gen.init_fused(size);
    EXPECT_EQ(gen.node_infos.size(), 2*4);
    optimizedSimulation sim(&boundaries,&gen);
    sim.run(0);
}
/**
 * Tests weather or not we can create a staircase with the functions.
 * @test
 * @see boundaryPointConstructor::steps_direction()
 */
TEST(InitTests, staircase_11) {
    unsigned int size = 8;
    point_t c = {size,size};
    point_t z = {0,0};
    boundaryPointConstructor boundaries(c);
    // boundaries.init_sliding_lid_side_chopped({20,10},30);
    boundaries.init_structure();
    boundaries.steps_direction(7,{1,1},&z,NO_BOUNDARY);
    // boundaries.visualize_2D_boundary(8);
    // in the creation logic the last point is not set so -1
    EXPECT_EQ(boundaries.boundary_structures[0]->boundary_points.size(), 8 +7 -1);
}

/**
 * Test weather or not a specific test boundary point constuction is ok.
 * @test
 */
TEST(InitTests, up_down_boundary) {
    unsigned int size = 10;
    point_t c = {size,size};
    point_t e1 = {1,0};
    point_t e2 = {4,4};
    point_t setter = {0,0};
    boundaryPointConstructor boundaries(c);
    boundaries.init_structure();
    boundaries.one_direction(8,{0,1},&setter,BOUNCE_BACK);
    boundaries.steps_direction(3,{1,1},&e1,BOUNCE_BACK);
    setter = {4,3};
    boundaries.set_point(&setter,BOUNCE_BACK);
    boundaries.steps_direction(3,{-1,1},&e2,BOUNCE_BACK);
    setter = {1,7};
    boundaries.set_point(&setter,BOUNCE_BACK);
    // boundaries.visualize_2D_boundary();
    EXPECT_EQ(boundaries.total_boundary_nodes(),22);
}

/**
 * Test/showcase of the surface variant of the boundary construction used for IBM_OUTER boundaries.
 * @test
 */
TEST(InitTests, basic_surface_test) {
    straight_t input;
    long canvas_size = 10;
    straightGenerator sg;
    // we put in a quader
    input.point = {1,1};
    input.direction = {0,1};
    input.max_t = 5;
    sg.add_surface(input);
    input.point = {1,6};
    input.direction = {1,0};
    input.max_t = 5;
    sg.add_surface(input);
    input.point = {6,6};
    input.direction = {0,-1};
    input.max_t = 5;
    sg.add_surface(input);
    input.point = {6,1};
    input.direction = {-1,0};
    input.max_t = 5;
    sg.add_surface(input);
    sg.surface_mass_center();
    EXPECT_EQ(sg.surfaces.size(),4);
    // node generator stuff
    nodeGenerator ng(&sg);
    ng.init_surface(canvas_size,0);
    // ng.visualize_2D_nodes();
    EXPECT_EQ(ng.node_infos.size(),16);
}

/**
 * Similar to basic_surface_test but the surface got moved.
 * @note this can be seen with a more tighly fit basic_surface_test
 * @test
 */
TEST(InitTests, basic_moved_surface) {
    point_t starter = {2.4,2.4};
    straight_t input;
    straightGenerator sg;
    long canvas_size = 10;
    double side_length = 4; // with a distance of 0.75 we should get 80 markers
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
    sg.surface_mass_center();
    EXPECT_EQ(sg.surfaces.size(),4);
    // node generator stuff
    nodeGenerator ng(&sg);
    ng.init_surface(canvas_size,1);
    // ng.visualize_2D_nodes();
    EXPECT_EQ(ng.node_infos.size(),16 + 20);
}

/**
 * We check on site ibm flagging.
 * @test
 */
TEST(InitTests, ibm_flagging) {
    point_t starter = {3.1,3.1};
    straight_t input;
    straightGenerator sg;
    long canvas_size = 15;
    double side_length = 6; // with a distance of 0.75 we should get 80 markers
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
    sg.surface_mass_center();
    EXPECT_EQ(sg.surfaces.size(),4);
    // node generator stuff
    nodeGenerator ng(&sg);
    ng.init_surface(canvas_size,1);
    // ng.visualize_2D_nodes();
    // check amount of ibm flagged nodes (not really necessary they can also be flagged as no boundary the markers hold all the boundary info
    int ibm_count = 0;
    for(auto node : ng.node_infos) {
        if(node->boundary == IBM_OUTER || node->boundary == IBM_INNER) {
            ++ibm_count;
        }
    }
    EXPECT_EQ(4*(5+7),ibm_count);
}

/**
 * We check that when we are on the surface we get an extetended marking.
 * @test
 */
TEST(InitTests, on_point_surface) {
    point_t starter = {3,3};
    straight_t input;
    straightGenerator sg;
    long canvas_size = 15;
    double side_length = 6; // with a distance of 0.75 we should get 80 markers
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
    sg.surface_mass_center();
    EXPECT_EQ(sg.surfaces.size(),4);
    // node generator stuff
    nodeGenerator ng(&sg);
    ng.init_surface(canvas_size,1);
    // ng.visualize_2D_nodes();
    // check amount of ibm flagged nodes (not really necessary they can also be flagged as no boundary the markers hold all the boundary info
    int ibm_count = 0;
    for(auto node : ng.node_infos) {
        if((node->boundary == IBM_OUTER) || (node->boundary == IBM_INNER)) {
            ++ibm_count;
        }
    }
    EXPECT_EQ(4*(4+6+8),ibm_count);
}

/*
// python stuff
def periodic_boundary_with_pressure_variations(grid,rho_in,rho_out):
    # get all the values
    rho, ux, uy = caluculate_real_values(grid)
    equilibrium = equilibrium_on_array(rho, ux, uy)
    ##########
    equilibrium_in = equilibrium_on_array(rho_in, ux[-2,:], uy[-2, :])
    # inlet 1,5,8
    grid[:, 0, :] = equilibrium_in + (grid[:, -2, :] - equilibrium[:, -2, :])

    # outlet 3,6,7
    equilibrium_out = equilibrium_on_array(rho_out, ux[1, :], uy[1, :])
    # check for correct sizes
    grid[:, -1, :] = equilibrium_out + (grid[:, 1, :] - equilibrium[:, 1, :])
*/