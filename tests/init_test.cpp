
#include "simulation.h"
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

TEST(InitTests,simulation_init_basic_test) {
    // init the boundary points -> then init the simulation ( and with it the node generator)
    // this is just a run test to see weather or not anything crashes durin 1 run
    int size = 5;
    point_t p = {size,size};
    boundaryPointConstructor boundaries(p);
    boundaries.init_quader();
    //
    simulation sim(&boundaries);
    sim.init();
    for(int i = 0; i< 1; ++i)
        sim.run();
    //sim.get_data(false);
}

TEST(InitTests, slidingLidBoundaries) {
    int size = 5;
    point_t p = {size,size};
    boundaryPointConstructor boundaries(p);
    boundaries.init_sliding_lid();
    double limit_y = boundaries.limits.y();
    for(auto b : boundaries.boundary_points) {
        if(b->point.y() == limit_y) {
            EXPECT_EQ(b->type, BOUNCE_BACK_MOVING);
        }
        else {
            EXPECT_EQ(b->type, BOUNCE_BACK);
        }
    }
}