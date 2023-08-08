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
#include "drawn_image_surface.h"
#include <gtest/gtest.h>

/**
 * Tests some of the functionality of the ranging key lookup.
 * @test
 */
TEST(FunctionalTest, ranging_key_look) {
    //
    double range = 2;
    // Fill a ranging point key hash
    std::vector<point_t> points;
    handle_t runner = 0;
    rangingPointKeyHash rpkh;
    for(int x = 0; x < 5; ++x) {
        for(int y = 0; y < 5; ++y) {
            point_t temp = {x,y};
            points.push_back(temp);
            rpkh.fill_key(runner,temp);
            runner++;
        }
    }
    // do a ranging scan
    point_t middle  = {2,2};
    std::vector<handle_t> back = rpkh.ranging_key_translation(middle, range);
    EXPECT_EQ(back.size(),points.size());
    EXPECT_TRUE(rpkh.ranging_key_look_for_specific(middle,range,12));
    // add an additional point
    point_t additional = {3.2,3.2};
    rpkh.fill_key(runner,additional);
    EXPECT_TRUE(rpkh.ranging_key_look_for_specific(middle,2,runner));
}

/**
 * We generate markers from a surface generated from an image, and put them into a neighborhood search thingi.
 * @note i noticed that wider placing helps with bug hunting, not usable for actual ibm calcs !
 * @test
 */
TEST(FunctionalTest, markers_on_surface_too_wide) {
    // node generator variables
    long canvas_size = 50;
    double marker_distance = 2; // intentionally to wide for actual ibm but easier to track
    kernelType_t kernel = KERNEL_C;
    double ibm_distance = kernel_id_to_lattice_search(kernel);
    // Load the image
    auto test_image = get_base_path();
    test_image.append("tests");
    test_image.append("test_images");
    test_image.append("cub.png");
    // call the drawer
    surfaceDrawer s(test_image);
    // std::vector<int> sel = {0,4,8,12};
    std::vector<int> sel = {0};
    s.run_selective(sel);
    s.surface_storage.surface_mass_center();
    nodeGenerator ng(&s.surface_storage);
    ng.init_surface_return(canvas_size,ibm_distance,marker_distance);
    // as a hard test prob best to write the markers into an
    rangingPointKeyHash rpkh;
    // loop over the marker points and fill the neigborhood search thingi
    handle_t h = 0;
    for(auto m : ng.markers->marker_points) {
        rpkh.fill_key(h,*m);
    }
    // find the ones and check distances
}
