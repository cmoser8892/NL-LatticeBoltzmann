// Include a library file to make sure proper includes are set
#include "one_step_simulation.h"
#include "node.h"
#include "functions.h"
#include "helper_functions.h"
#include "neighborhood.h"
#include "image_converter.h"
#include <gtest/gtest.h>

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
    // testing the equilbirum function
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
    // exceptions
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
    optimizedSimulation sim(&boundaries, &gen);
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
    sim.offset_node = (step & 0x1) * 9;
    for(auto n : sim.nodes) {
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

TEST(FunctionalTest, simple_use_case_straight) {
    // simple quader without any special cases for the straight surface generator
    int step = 0;
    int size = 4;
    point_t start = {0,0};
    point_t end = {size,0};
    point_t sim_area = {size,size};
    // init
    boundaryPointConstructor boundaries(sim_area);
    boundaries.init_quader();
    // boundaries.visualize_2D_boundary();
    // straight tests
    straightGenerator s(&boundaries);
    s.init();
    EXPECT_EQ(s.surfaces.size(),4);
    for(auto surface : s.surfaces) {
        EXPECT_EQ(surface->max_t, 3);
    }
}

TEST(FunctionalTest, straight_inner_missing_57) {
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
    //
    straightGenerator s(&boundaries);
    s.init();
    // check the surfaces
    EXPECT_EQ(s.surfaces.size(),8);
    // control the values
    int ones = 0;
    int tows = 0;
    int errors = 0;
    for(auto surface : s.surfaces) {
        if(surface->max_t == 1) {
            ++ones;
        }
        else if(surface->max_t == 2) {
            ++tows;
        }
        else {
            ++errors;
        }
    }
    EXPECT_EQ(ones,4);
    EXPECT_EQ(tows,4);
    EXPECT_EQ(errors,0);
}

TEST(FunctionalTest, straight_inner_missing_68) {
    // sizes matters
    int step = 0;
    int size = 4;
    point_t start = {0,0};
    point_t end = {size,0};
    point_t sim_area = {size,size};
    boundaryPointConstructor boundaries(sim_area);
    // init boundaries
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
    //
    straightGenerator s(&boundaries);
    s.init();
    // check the surfaces
    EXPECT_EQ(s.surfaces.size(),8);
    // control the values
    int ones = 0;
    int tows = 0;
    int errors = 0;
    for(auto surface : s.surfaces) {
        if(surface->max_t == 1) {
            ++ones;
        }
        else if(surface->max_t == 2) {
            ++tows;
        }
        else {
            ++errors;
        }
    }
    EXPECT_EQ(ones,4);
    EXPECT_EQ(tows,4);
    EXPECT_EQ(errors,0);
}


TEST(FunctionalTest, straight_deletions) {
    unsigned int size = 10;
    point_t c = {size,size};
    point_t e1 = {1,0};
    point_t e2 = {4,4};
    point_t setter = {0,0};
    // set up a boundary test structure to tests interrupting lines
    boundaryPointConstructor boundaries(c);
    boundaries.init_structure();
    boundaries.one_direction(8,{0,1},&setter,BOUNCE_BACK);
    boundaries.steps_direction(3,{1,1},&e1,BOUNCE_BACK);
    setter = {4,3};
    boundaries.set_point(&setter,BOUNCE_BACK);
    boundaries.steps_direction(3,{-1,1},&e2,BOUNCE_BACK);
    setter = {1,7};
    boundaries.set_point(&setter,BOUNCE_BACK);
    EXPECT_EQ(boundaries.total_boundary_nodes(),22);
    // boundaries.visualize_2D_boundary();
    // manually set up the straight generator to tests out performance
    straightGenerator s(&boundaries);
    s.init();
    EXPECT_EQ(s.surfaces.size(),12);
    // control the values
    int ones = 0;
    int twos = 0;
    int threes = 0;
    int sevens = 0;
    int errors = 0;
    for(auto surface : s.surfaces) {
        if(surface->max_t == 1) {
            ++ones;
        }
        else if(surface->max_t == 2) {
            ++twos;
        }
        else if(surface->max_t == 3) {
            ++threes;
        }
        else if(surface->max_t == 7) {
            ++sevens;
        }
        else {
            ++errors;
        }
    }
    EXPECT_EQ(ones,8);
    EXPECT_EQ(twos,2);
    EXPECT_EQ(threes, 1);
    EXPECT_EQ(sevens,1);
    EXPECT_EQ(errors,0);
    // test out the fully integrated thing in the node generator
    nodeGenerator gen(&boundaries);
    gen.init(size);
    // gen.visualize_2D_nodes();
    EXPECT_EQ(gen.node_infos.size(), 22 + 6 + 4 + 2);
}

TEST(FunctionalTest, multiple_interruptions) {
    unsigned int size = 10;
    point_t c = {4,size};
    point_t setter = {0,0};
    // set up a boundary test structure to tests interrupting lines
    boundaryPointConstructor boundaries(c);
    boundaries.init_structure();
    boundaries.one_direction(10,{0,1},&setter,BOUNCE_BACK);
    setter = {1,0};
    boundaries.steps_direction(2,{1,1},&setter,BOUNCE_BACK);
    setter = {3,2};
    boundaries.steps_direction(2,{-1,1},&setter,BOUNCE_BACK);
    setter = {2,5};
    boundaries.steps_direction(1,{1,1},&setter,BOUNCE_BACK);
    setter = {3,6};
    boundaries.steps_direction(2,{-1,1},&setter,BOUNCE_BACK);
    setter = {2,9};
    boundaries.steps_direction(1,{-1,-1},&setter,BOUNCE_BACK);
    // visualizer
    // boundaries.visualize_2D_boundary();
    straightGenerator s(&boundaries);
    s.init();
    // control
    int ones = 0;
    int twos = 0;
    int nines = 0;
    int errors = 0;
    for(auto surface : s.surfaces) {
        if(surface->max_t == 1) {
            ++ones;
        }
        else if(surface->max_t == 2) {
            ++twos;
        }
        else if(surface->max_t == 9) {
            ++nines;
        }
        else {
            ++errors;
        }
    }
    EXPECT_EQ(ones,5);
    EXPECT_EQ(twos,6);
    EXPECT_EQ(nines,1);
    EXPECT_EQ(errors,0);
    EXPECT_EQ(s.surfaces.size(),12);
    // node generator master test
    nodeGenerator gen(&boundaries);
    gen.init(size);
    // gen.visualize_2D_nodes();
    EXPECT_EQ(gen.node_infos.size(), boundaries.total_boundary_nodes() + 8 + 2);
}

TEST(FunctionalTest,comparision) {
    point_t t = {2,1};
    EXPECT_TRUE(compare_two_points(&t,&t));
}

TEST(FunctionalTest, straight_bump) {
    unsigned int size = 6;
    point_t c = {size,size};
    point_t setter = {2,1};
    boundaryPointConstructor boundaries(c);
    boundaries.init_quader();
    boundaries.set_point(&setter,BOUNCE_BACK);
    setter = {3,1};
    boundaries.set_point(&setter,BOUNCE_BACK);
    setter = {2,2};
    boundaries.set_point(&setter,BOUNCE_BACK);
    setter = {3,2};
    boundaries.set_point(&setter,BOUNCE_BACK);
    // boundaries.visualize_2D_boundary();
    // straight test
    straightGenerator s(&boundaries);
    s.init();
    EXPECT_EQ(s.surfaces.size(),8);
    // control the values
    int ones = 0;
    int twos = 0;
    int fives = 0;
    int errors = 0;
    for(auto surface : s.surfaces) {
        if(surface->max_t == 1) {
            EXPECT_EQ(surface->point.x(),2);
            EXPECT_EQ(surface->point.y(),2);
            ++ones;
        }
        else if(surface->max_t == 2) {
            ++twos;
        }
        else if(surface->max_t == 5) {
            ++fives;
        }
        else {
            ++errors;
        }
    }
    EXPECT_EQ(ones,1);
    EXPECT_EQ(twos,4);
    EXPECT_EQ(fives, 3);
    EXPECT_EQ(errors,0);
    // direct node gen test
    nodeGenerator gen(&boundaries);
    gen.init_fused(size);
    // gen.visualize_2D_nodes();
    EXPECT_EQ(gen.node_infos.size(),4 +8);
}

TEST(FunctionalTest, multiple_bumps) {
    unsigned int size = 9;
    point_t c = {size,size};
    boundaryPointConstructor boundaries(c);
    boundaries.init_quader();
    point_t setter = {2,1};
    boundaries.set_point(&setter,BOUNCE_BACK);
    setter = {3,1};
    boundaries.set_point(&setter,BOUNCE_BACK);
    setter = {2,2};
    boundaries.set_point(&setter,BOUNCE_BACK);
    setter = {3,2};
    boundaries.set_point(&setter,BOUNCE_BACK);
    setter = {5,1};
    boundaries.set_point(&setter,BOUNCE_BACK);
    setter = {6,1};
    boundaries.set_point(&setter,BOUNCE_BACK);
    setter = {5,2};
    boundaries.set_point(&setter,BOUNCE_BACK);
    setter = {6,2};
    boundaries.set_point(&setter,BOUNCE_BACK);
    setter = {2,7};
    boundaries.set_point(&setter,BOUNCE_BACK);
    setter = {3,7};
    boundaries.set_point(&setter,BOUNCE_BACK);
    setter = {2,6};
    boundaries.set_point(&setter,BOUNCE_BACK);
    setter = {3,6};
    boundaries.set_point(&setter,BOUNCE_BACK);
    // boundaries.visualize_2D_boundary();
    nodeGenerator gen(&boundaries);
    gen.init_fused(size);
    // gen.visualize_2D_nodes();
    EXPECT_EQ(gen.node_infos.size(),7+3+3+7+5+5+7);
}

TEST(FunctionalTest, special_case_little_bump) {
    unsigned int size = 20;
    unsigned int sub_size = 9;
    point_t c = {size,size};
    point_t s = {sub_size,sub_size};
    boundaryPointConstructor boundaries(c);
    vector_t shift = {4,4};
    boundaries.init_quader(shift,s);
    point_t setter = {2,1};
    setter += shift;
    boundaries.set_point(&setter,BOUNCE_BACK);
    setter = {3,1};
    setter += shift;
    boundaries.set_point(&setter,BOUNCE_BACK);
    setter = {2,2};
    setter += shift;
    boundaries.set_point(&setter,BOUNCE_BACK);
    setter = {3,2};
    setter += shift;
    boundaries.set_point(&setter,BOUNCE_BACK);
    setter = {5,1};
    setter += shift;
    boundaries.set_point(&setter,BOUNCE_BACK);
    setter = {6,1};
    setter += shift;
    boundaries.set_point(&setter,BOUNCE_BACK);
    setter = {5,2};
    setter += shift;
    boundaries.set_point(&setter,BOUNCE_BACK);
    setter = {6,2};
    setter += shift;
    boundaries.set_point(&setter,BOUNCE_BACK);
    setter = {2,7};
    setter += shift;
    boundaries.set_point(&setter,BOUNCE_BACK);
    setter = {3,7};
    setter += shift;
    boundaries.set_point(&setter,BOUNCE_BACK);
    setter = {2,6};
    setter += shift;
    boundaries.set_point(&setter,BOUNCE_BACK);
    setter = {3,6};
    setter += shift;
    boundaries.set_point(&setter,BOUNCE_BACK);
    // boundaries.visualize_2D_boundary();
    nodeGenerator gen(&boundaries);
    gen.init_fused(size);
    // gen.visualize_2D_nodes();
    EXPECT_EQ(gen.node_infos.size(),7+3+3+7+5+5+7);
}

TEST(FunctionalTest, one_high_side_walls) {
    int size = 6;
    point_t sim_area = {size,size};
    // init
    boundaryPointConstructor boundaries(sim_area);
    boundaries.init_quader();
    point_t setter = {2,1};
    boundaries.one_direction(2,{1,0},&setter,BOUNCE_BACK);
    // boundaries.visualize_2D_boundary();
    nodeGenerator gen(&boundaries);
    gen.init_fused(size);
    // gen.visualize_2D_nodes();
    EXPECT_EQ(gen.node_infos.size(),4+3+3+4);
}

TEST(FunctionalTest, wierd_bump_top) {
    // init variables
    unsigned int size = 20;
    unsigned int sub_size = 9;
    point_t c = {size,size};
    point_t s = {sub_size,sub_size};
    vector_t shift = {4,4};
    point_t setter = {2,1};
    // basic setup
    boundaryPointConstructor boundaries(c);
    boundaries.init_structure();
    // structural init
    setter = {4,4};
    boundaries.one_direction(10,{0,1},&setter,BOUNCE_BACK);
    setter = {4,14};
    boundaries.one_direction(10,{1,0},&setter,BOUNCE_BACK);
    setter = {14,14};
    boundaries.one_direction(10,{0,-1},&setter,BOUNCE_BACK);
    // delicate stuff
    setter = {5,4};
    boundaries.one_direction(3,{1,0},&setter,BOUNCE_BACK);
    setter = {7,5};
    boundaries.one_direction(3,{0,1},&setter,BOUNCE_BACK);
    setter = {7,8};
    boundaries.one_direction(3,{1,0},&setter,BOUNCE_BACK);
    setter = {7,9};
    boundaries.one_direction(2,{1,0},&setter,BOUNCE_BACK);  /// defining
    setter = {10,8};
    boundaries.one_direction(4,{0,-1},&setter,BOUNCE_BACK);
    setter = {10,4};
    boundaries.one_direction(5,{1,0},&setter,BOUNCE_BACK);
    // boundaries.visualize_2D_boundary();
    // it is always a good idea to look at straight.temporary for a good test
    // test for right surface
    straightGenerator stg(&boundaries);
    stg.init();
    EXPECT_EQ(stg.surfaces.size(),10);
    // testing surface lengths + individually inspect
    // control the values
    int ones = 0;
    int twos = 0;
    int threes = 0;
    int fours = 0;
    int fives = 0;
    int tens = 0;
    int errors = 0;
    for(auto surface : stg.surfaces) {
        switch(int(surface->max_t)) {
        case 1:
            ++ones;
            break;
        case 2:
            ++twos;
            break;
        case 3:
            ++threes;
            break;
        case 4:
            ++fours;
            break;
        case 5:
            ++fives;
            break;
        case 10:
            ++tens;
            break;
        default:
            ++errors;
            break;
        }

    }
    EXPECT_EQ(ones,2);
    EXPECT_EQ(twos,1);
    EXPECT_EQ(threes,1);
    EXPECT_EQ(fours,2);
    EXPECT_EQ(fives,1);
    EXPECT_EQ(tens,3);
    EXPECT_EQ(errors,0);
}

TEST(FunctionalTest, wierd_bump_bottom) {
    // init variables
    unsigned int size = 20;
    unsigned int sub_size = 9;
    point_t c = {size,size};
    point_t s = {sub_size,sub_size};
    vector_t shift = {4,4};
    point_t setter = {2,1};
    // basic setup
    boundaryPointConstructor boundaries(c);
    boundaries.init_structure();
    // structural init
    setter = {4,4};
    boundaries.one_direction(10,{0,1},&setter,BOUNCE_BACK);
    setter = {4,14};
    boundaries.one_direction(10,{1,0},&setter,BOUNCE_BACK);
    setter = {14,14};
    boundaries.one_direction(10,{0,-1},&setter,BOUNCE_BACK);
    // delicate stuff
    setter = {5,4};
    boundaries.one_direction(3,{1,0},&setter,BOUNCE_BACK);
    setter = {7,5};
    boundaries.one_direction(3,{0,1},&setter,BOUNCE_BACK);
    setter = {7,8};
    boundaries.one_direction(3,{1,0},&setter,BOUNCE_BACK);
    setter = {9,9};
    boundaries.one_direction(2,{1,0},&setter,BOUNCE_BACK);  /// defining
    setter = {10,8};
    boundaries.one_direction(4,{0,-1},&setter,BOUNCE_BACK);
    setter = {10,4};
    boundaries.one_direction(5,{1,0},&setter,BOUNCE_BACK);
    // boundaries.visualize_2D_boundary();
    // it is always a good idea to look at straight.temporary for a good test
    // test for right surface
    straightGenerator stg(&boundaries);
    stg.init();
    EXPECT_EQ(stg.surfaces.size(),10);
    // testing surface lengths + individually inspect
    // control the values
    int ones = 0;
    int twos = 0;
    int threes = 0;
    int fours = 0;
    int fives = 0;
    int tens = 0;
    int errors = 0;
    for(auto surface : stg.surfaces) {
        switch(int(surface->max_t)) {
        case 1:
            ++ones;
            break;
        case 2:
            ++twos;
            break;
        case 3:
            ++threes;
            break;
        case 4:
            ++fours;
            break;
        case 5:
            ++fives;
            break;
        case 10:
            ++tens;
            break;
        default:
            ++errors;
            break;
        }

    }
    EXPECT_EQ(ones,2);
    EXPECT_EQ(twos,1);
    EXPECT_EQ(threes,1);
    EXPECT_EQ(fours,2);
    EXPECT_EQ(fives,1);
    EXPECT_EQ(tens,3);
    EXPECT_EQ(errors,0);
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
    // ic.boundaries->visualize_2D_boundary();
    nodeGenerator gen(ic.boundaries);
    // if the fused init runs this test is considered complete
    gen.init_fused(ic.return_basic_size());
    // gen.visualize_2D_nodes();
    // nessessary to to a visual check
    EXPECT_TRUE(true);
}

TEST(FunctionalTest,InnerCorner) {
    // make sure that sth like that doesnt crash the node generator
    // init variables
    unsigned int size = 20;
    unsigned int sub_size = 9;
    point_t c = {size,size};
    point_t s = {sub_size,sub_size};
    vector_t shift = {4,4};
    point_t setter = {2,1};
    // basic setup
    boundaryPointConstructor boundaries(c);
    boundaries.init_structure();
    // structural init
    setter = {4,4};
    boundaries.one_direction(10,{0,1},&setter,BOUNCE_BACK);
    setter = {4,14};
    boundaries.one_direction(10,{1,0},&setter,BOUNCE_BACK);
    setter = {14,14};
    boundaries.one_direction(10,{0,-1},&setter,BOUNCE_BACK);
    // delicate stuff
    setter = {5,4};
    boundaries.one_direction(3,{1,0},&setter,BOUNCE_BACK);
    setter = {7,5};
    boundaries.one_direction(3,{0,1},&setter,BOUNCE_BACK);
    setter = {7,8};
    boundaries.one_direction(3,{1,0},&setter,BOUNCE_BACK);
    setter = {7,9};
    boundaries.one_direction(2,{1,0},&setter,BOUNCE_BACK);  /// defining
    // setter = {10,9};
    // boundaries.one_direction(2,{1,0},&setter,BOUNCE_BACK);  /// defining
    setter = {10,8};
    boundaries.one_direction(4,{0,-1},&setter,BOUNCE_BACK);
    setter = {10,4};
    boundaries.one_direction(5,{1,0},&setter,BOUNCE_BACK);
    // boundaries.visualize_2D_boundary();
    nodeGenerator gen(&boundaries);
    gen.init_fused(size);
    // gen.visualize_2D_nodes();
    // as long as the test does not crash it is considered ok
    EXPECT_TRUE(true);
}

TEST(FunctionalTest, StarIntersections) {
    unsigned int size = 20;
    unsigned int sub_size = 9;
    point_t c = {size,size};
    point_t s = {sub_size,sub_size};
    vector_t shift = {4,4};
    point_t setter = {2,1};
    // basic setup
    boundaryPointConstructor boundaries(c);
    boundaries.init_structure();
    // structural init
    setter = {4,4};
    boundaries.one_direction(10,{0,1},&setter,BOUNCE_BACK);
    setter = {4,14};
    boundaries.one_direction(10,{1,0},&setter,BOUNCE_BACK);
    setter = {14,14};
    boundaries.one_direction(10,{0,-1},&setter,BOUNCE_BACK);
    // delicate stuff
    setter = {5,4};
    boundaries.one_direction(3,{1,0},&setter,BOUNCE_BACK);
    setter = {7,5};
    boundaries.one_direction(3,{0,1},&setter,BOUNCE_BACK);
    setter = {7,8};
    boundaries.one_direction(3,{1,0},&setter,BOUNCE_BACK);
    setter = {7,9};
    boundaries.one_direction(2,{1,0},&setter,BOUNCE_BACK);  /// defining
    // setter = {10,9};
    // boundaries.one_direction(2,{1,0},&setter,BOUNCE_BACK);  /// defining
    setter = {10,8};
    boundaries.one_direction(4,{0,-1},&setter,BOUNCE_BACK);
    setter = {10,4};
    boundaries.one_direction(5,{1,0},&setter,BOUNCE_BACK);
    // boundaries.visualize_2D_boundary();
    nodeGenerator gen(&boundaries);
    gen.init(size);
    // gen.visualize_2D_nodes();
}
