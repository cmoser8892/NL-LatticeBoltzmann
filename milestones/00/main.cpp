
/*
 * vallgrind call
valgrind --tool=callgrind --dump-instr=yes (p)
 */

#include "helper_functions.h"
#include "drawn_image_surface.h"

int main() {
    // Load the image
    auto test_image = get_base_path();
    test_image.append("tests");
    test_image.append("test_images");
    test_image.append("donut.png");
    // call the drawer
    surfaceDrawer s(test_image);
    s.run();
    // end
    return 0;
}
