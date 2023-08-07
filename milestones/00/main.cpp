
/*
 * vallgrind call
valgrind --tool=callgrind --dump-instr=yes (p)
 */

#include "helper_functions.h"
#include "drawn_image_surface.h"

/**
 * 00 main, testing stuff.
 * @note good rule of thumb is all 4 contours need to be selected
 * @return
 */
int main() {
    // Load the image
    auto test_image = get_base_path();
    test_image.append("tests");
    test_image.append("test_images");
    test_image.append("blobs.png");
    // call the drawer
    surfaceDrawer s(test_image);
    std::vector<int> sel = {0,4,8,12};
    s.run_selective(sel);
    // end
    return 0;
}
