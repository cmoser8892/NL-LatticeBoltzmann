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
    EXPECT_EQ(boundaries.total_boundary_nodes(), (size-1)*(size-1));
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
    for(auto bs : boundaries.boundary_structures) {
        for( auto b: bs->boundary_points) {
            if(b->point.y() == size -1) {
                EXPECT_EQ(b->type, BOUNCE_BACK_MOVING);
            }
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
    EXPECT_EQ(boundaries.total_boundary_nodes(),4*(sub_size-1));
    EXPECT_EQ(st.surfaces.size(),4*(sub_size-1));
    EXPECT_EQ(st.surfaces.size(),boundaries.total_boundary_nodes());
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
    EXPECT_EQ(boundaries.total_boundary_nodes(), (sub_size-1)*4);
    // check if even with the bulge the sizes are still the same
    EXPECT_EQ(boundaries.total_boundary_nodes(),(sub_size-1)*4);
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
    EXPECT_EQ(boundaries.total_boundary_nodes(), ((outer_size-1) + (inner_size-1))*4);
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
    // still kinda spooked will try out the write out of the rho field
    unsigned int size = 30;
    unsigned int sub_size = 20;
    point_t c = {size,size};
    point_t p = {sub_size +5,sub_size};
    point_t k = {4,10};
    boundaryPointConstructor boundaries(p);
    // boundaries.init_sliding_lid_side_chopped({20,10},30);
    boundaries.init_sliding_lid_inner({3,5},{9,7},k);
    // boundaries.visualize_2D_boundary(30);
    nodeGenerator gen(&boundaries);
    gen.init(size);
    // gen.visualize_2D_nodes(30);
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

TEST(InitTests, init_out_inner_rho_writeout) {
    // checks number of rho zeros
    unsigned int size = 30;
    unsigned int sub_size = 20;
    point_t c = {size,size};
    point_t p = {sub_size,sub_size+2};
    point_t k = {5,6};
    boundaryPointConstructor boundaries(p);
    // boundaries.init_sliding_lid_side_chopped({20,10},30);
    boundaries.init_sliding_lid_inner({3,5},{9,7},k);
    nodeGenerator gen(&boundaries);
    gen.init(size);
    simulation sim(&boundaries,&gen);
    sim.init();
    flowfield_t rho;
    rho.resize(size,size);
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

TEST(InitTests, inner_outer_master_test) {
    unsigned int size = 10;
    unsigned int sub_size = 8;
    unsigned int white_list = 6;
    unsigned int inner_size = 4;
    point_t c = {size,size}; // size or the canvas
    point_t p = {sub_size,sub_size}; // size of the  quader
    point_t k = {inner_size,inner_size}; // size of the inner one
    boundaryPointConstructor boundaries(p);
    // boundaries.init_sliding_lid_side_chopped({20,10},30);
    boundaries.init_sliding_lid_inner({1,1},{3,3},k);
    // boundaries.visualize_2D_boundary(size);
    EXPECT_EQ(boundaries.total_boundary_nodes(),(sub_size-1)*4 + (inner_size -1)*4);
    nodeGenerator gen(&boundaries);
    // check boundaries right at least
    gen.init(size);
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

TEST(InitTests, fused_inner_outer_init_simple) {
    unsigned int size = 10;
    unsigned int sub_size = 8;
    unsigned int white_list = 6;
    unsigned int inner_size = 4;
    point_t c = {size,size}; // size or the canvas
    point_t p = {sub_size,sub_size}; // size of the  quader
    point_t k = {inner_size,inner_size}; // size of the inner one
    boundaryPointConstructor boundaries(p);
    // boundaries.init_sliding_lid_side_chopped({20,10},30);
    boundaries.init_sliding_lid_inner({1,1},{3,3},k);
    // boundaries.visualize_2D_boundary(size);
    EXPECT_EQ(boundaries.total_boundary_nodes(),(sub_size-1)*4 + (inner_size -1)*4);
    nodeGenerator gen(&boundaries);
    // check boundaries right at least
    gen.init_fused(size);
    // gen.visualize_2D_nodes(size);
    int expected_total_node_number = (p.x()-2)*(p.y()-2) - (k.x()*k.y());
    EXPECT_EQ(gen.node_infos.size(),expected_total_node_number);
}

TEST(InitTests, fused_inner_outer_init_shifted) {
    // checks number of rho zeros
    unsigned int size = 30;
    unsigned int sub_size = 20;
    point_t c = {size,size};
    point_t p = {sub_size,sub_size+2};
    point_t k = {5,6};
    boundaryPointConstructor boundaries(p);
    // boundaries.init_sliding_lid_side_chopped({20,10},30);
    boundaries.init_sliding_lid_inner({3,5},{9,7},k);
    nodeGenerator gen(&boundaries);
    gen.init_fused(size);
    simulation sim(&boundaries,&gen);
    sim.fused_init();
    // check
    int expected_total_node_number = (p.x()-2)*(p.y()-2) - (k.x()*k.y());
    EXPECT_EQ(sim.nodes.size(),expected_total_node_number);
}

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
    EXPECT_EQ(pb, 20+20);
}