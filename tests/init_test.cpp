
#include "simulation.h"
#include "helper_functions.h"
#include "functions.h"
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


TEST(InitTests,simulation_init_run) {
    // init the boundary points -> then init the simulation ( and with it the node generator)
    // this is just a run test to see weather or not anything crashes durin 1 run
    int size = 3;
    int steps = 5;
    point_t p = {size,size};
    boundaryPointConstructor boundaries(p);
    boundaries.init_quader();
    //
    simulation sim(&boundaries);
    sim.init();
    debug_node(sim.nodes.at(0),true);
    for(int i = 0; i < steps; ++i) {
        sim.run();
        debug_node(sim.nodes.at(0),true);
    }
    sim.get_data(false);
}