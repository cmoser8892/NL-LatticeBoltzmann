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
 * Tests the readin in the new drawer class.
 * @test
 */
TEST(FunctionalTest,drawer_basics) {
    // test image setupa
    auto test_image = get_base_path();
    test_image.append("tests");
    test_image.append("test_images");
    test_image.append("20_blocks.bmp");
    surfaceDrawer s(test_image);
    s.init();
    EXPECT_EQ(s.points.size(),20);
}
