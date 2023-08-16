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
     if(0) {
         std::cout << "Point " << s.point.x() <<" ," << s.point.y()
                   << "\n Direction:" << s.direction.x() << " ," << s.direction.y() << " Length: " << s.max_t << std::endl;
     }
 }

 /**
  * Loops throught the selector vector and returns true if the current contour is in there
  * @note selectors start at 0 i guess.
  * @param sel
  * @param current_contour
  * @return
  */
bool surfaceDrawer::selector(std::vector<int> sel, int current_contour) {
     bool returns = false;
     for(auto i : sel) {
        if(i == current_contour) {
            returns = true;
            break;
        }
     }
     return returns;
}

/**
 * Constructor.
 * @param p
 */
surfaceDrawer::surfaceDrawer(std::filesystem::path p) {
     path = p;
}

/**
 * Runs an opencv instance, only for the outer contour thou.
 * @details https://docs.opencv.org/4.x/d9/d8b/tutorial_py_contours_hierarchy.html
 * @note we just want the contours given in our image
 * @attention prob really specific to the given use case and how we draw in general.
 */
void surfaceDrawer::run() {
     std::vector<int> sel = {0};
     run_selective(sel);
}

/**
 * Runs an opencv instance, we can chose form the contours found which one we want in our final surface.
 * @note Trial and error is king here.
 * @param s
 */
void surfaceDrawer::run_selective(std::vector<int> s) {
     run_non_connecting(s,true);
}

/**
 *
 * @param s
 * @param connect_start_end
 */
void surfaceDrawer::run_non_connecting(std::vector<int> s, bool connect_start_end) {
     // read the image
     cv::Mat image = cv::imread(path.string(), cv::IMREAD_GRAYSCALE);
     // throw an runtime error in case we could not read that
     if (image.empty()) {
         throw std::runtime_error("Error: Unable to load the image." );
     }
     // detect the edges
     cv::Mat edges;
     cv::Canny(image, edges, 50, 150);
     // cv::imshow("Edges", edges);
     // cv::waitKey(0);
     // detect the contours
     std::vector<std::vector<cv::Point>> contours;
     // RETR_EXTERNAL vs Tree
     cv::findContours(edges, contours, cv::RETR_TREE, cv::CHAIN_APPROX_SIMPLE);
     // input into the straigth generator
     if (!contours.empty()) {
         // we only look at the top level surface
         // const std::vector<cv::Point> &contour = contours[0];
         // contour loop
         int i = -1;
         for(auto contour : contours) {
             ++i;
             if(selector(s,i)) {
                 size_t run = contour.size();
                 if(!connect_start_end) {
                     run--;
                 }
                 for (size_t i = 0; i < run; i++) {
                     int next_index = (i + 1) % contour.size();
                     // point translation
                     cv::Point first = contour[i];
                     cv::Point second = contour[next_index];
                     point_t intern_first = {first.x,first.y};
                     point_t intern_second = {second.x,second.y};
                     // input into the straight generator
                     add_surface(intern_first,intern_second);
                 }
             }
         }
     }
     surface_storage.write_out_surface();
}

/**
 * Closes the surface if we called run_non_connecting.
 * @param draw_size
 */
void surfaceDrawer::close_open_surface(vector_t draw_size) {
     std::vector<std::vector<point_t>> endpoints;
     endpoints.resize(4);
     for(auto s : surface_storage.surfaces) {
         // we need to check out where they go?
         for(int i = 0; i < s->point.size(); ++i) {
             point_t test_point = s->point + s->direction * s->max_t * i;
             uint8_t test_value = point_on_boarder(&test_point,&draw_size);
             if(test_value>0) {
                 endpoints[test_value-1].push_back(test_point);
             }
         }
     }
     for(auto b : endpoints) {
         if(!b.empty()) {
             if(b.size() >= 2) {
                 point_t start = b[0];
                 point_t end = b[b.size()-1];
                 vector_t direction = end -start;
                 straight_t bend;
                 bend.point = start;
                 bend.direction = direction.normalized();
                 bend.max_t = direction.norm();
                 surface_storage.add_surface(bend);
                 std::cout << "added" << std::endl;
             }
             else {
                 std::cerr << "Could not identify endpoints" << std::endl;
             }
         }
     }
}
