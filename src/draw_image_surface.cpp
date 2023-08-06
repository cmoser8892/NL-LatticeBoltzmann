#include <opencv2/opencv.hpp>
#include "drawn_image_surface.h"

/**
 * Adds a surface to the surface_storage, aka the straightGenerator.
 * @param current
 * @param previous
 */
 void surfaceDrawer::add_surface(point_t current, point_t previous) {
     vector_t direction = current - previous;
     straight_t  s;
     s.point = previous;
     s.direction = direction.normalized();
     s.max_t = direction.norm();
     surface_storage.add_surface(s);
     // printer
     if(1) {
         std::cout << "Point " << s.point.x() <<" ," << s.point.y()
                   << "\n Direction:" << s.direction.x() << " ," << s.direction.y() << " Length: " << s.max_t << std::endl;
     }
 }

/**
 * Constructor.
 * @param p
 */
surfaceDrawer::surfaceDrawer(std::filesystem::path p) {
     path = p;
}

/**
 * Runs an opencv instance.
 * @details https://docs.opencv.org/4.x/d9/d8b/tutorial_py_contours_hierarchy.html
 * @note we just want the contours given in our image
 * @attention prob really specific to the given use case and how we draw in general.
 */
void surfaceDrawer::run() {
     // read the image
     cv::Mat image = cv::imread(path.string(), cv::IMREAD_GRAYSCALE);
     // throw an runtime error in case we could not read that
     if (image.empty()) {
         throw std::runtime_error("Error: Unable to load the image." );
     }
     // detect the edges
     cv::Mat edges;
     cv::Canny(image, edges, 50, 150);
     // detect the contours
     std::vector<std::vector<cv::Point>> contours;
     cv::findContours(edges, contours, cv::RETR_EXTERNAL, cv::CHAIN_APPROX_SIMPLE);
     // input into the straigth generator
     if (!contours.empty()) {
         // we only look at the top level surface
         const std::vector<cv::Point> &contour = contours[0];
         // contour loop
         for (size_t i = 0; i < contour.size(); i++) {
             int nextIndex = (i + 1) % contour.size();
             // point translation
             cv::Point first = contour[i];
             cv::Point second = contour[nextIndex];
             point_t intern_first = {first.x,first.y};
             point_t intern_second = {second.x,second.y};
             if(i == 49) {
                 std::cout << "s" << std::endl;
             }
             // input into the straight generator
             add_surface(intern_first,intern_second);
         }
     }
     surface_storage.write_out_surface();
}
