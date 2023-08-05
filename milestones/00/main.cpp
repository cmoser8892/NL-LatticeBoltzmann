
/*
 * vallgrind call
valgrind --tool=callgrind --dump-instr=yes (p)
 */

#include <iostream>
#include <vector>
#include <opencv2/opencv.hpp>
#include "helper_functions.h"
#include "straight.h"

using namespace cv;
using namespace std;

// Function to apply the Douglas-Peucker line simplification
void simplifyLines(const vector<Point>& inputPoints, vector<Point>& outputPoints, double epsilon) {
    cout << "Input points: " << inputPoints.size() << endl;
    if (inputPoints.size() < 2) {
        outputPoints = inputPoints;
        return;
    }

    vector<Point> simplified;
    simplified.push_back(inputPoints.front());
    simplified.push_back(inputPoints.back());

    for (size_t i = 2; i < inputPoints.size(); i++) {
        double distance = pointPolygonTest(inputPoints, inputPoints[i], true);
        if (distance > epsilon) {
            simplified.push_back(inputPoints[i]);
        }
    }
    cout << "Simplified points: " << outputPoints.size() << endl;
    outputPoints = simplified;
}

int main() {
    // Load the image
    auto test_image = get_base_path();
    test_image.append("tests");
    test_image.append("test_images");
    test_image.append("extra_small.bmp");
    Mat image = imread(test_image.string(), IMREAD_GRAYSCALE);
    if (image.empty()) {
        cout << "Error: Unable to load the image." << endl;
        return -1;
    }
    // Apply edge detection
    Mat edges;
    Canny(image, edges, 50, 150);

    // Find contours from the edges
    vector<vector<Point>> contours;  // Ensure that Point data type is used here
    findContours(edges, contours, RETR_EXTERNAL, CHAIN_APPROX_SIMPLE);

    // Create a blank image to draw the surface
    Mat surfaceImage = Mat::zeros(image.size(), CV_8UC3);

    // Connect the points of the first contour to recreate the surface
    if (!contours.empty()) {
        const vector<Point>& contour = contours[0];
        for (size_t i = 0; i < contour.size(); i++) {
            int nextIndex = (i + 1) % contour.size();
            line(surfaceImage, contour[i], contour[nextIndex], Scalar(0, 0, 255), 2);
        }
    }

    // Display the regenerated surface
    imshow("Regenerated Surface", surfaceImage);
    waitKey(0);
    // own
    if(1) {
        // input into the straigth generator
        straightGenerator sg;
        if (!contours.empty()) {
            const vector<Point> &contour = contours[0];
            // contour loop
            for (size_t i = 0; i < contour.size(); i++) {
                int nextIndex = (i + 1) % contour.size();
                // point translation
                Point first = contour[i];
                Point second = contour[nextIndex];
                point_t intern_first = {first.x,first.y};
                point_t intern_second = {second.x,second.y};
                if(i == 49) {
                    std::cout << "s" << std::endl;
                }
                // input into the straight generator
                vector_t direction = intern_second - intern_first;
                straight_t s;
                s.point = intern_first;
                s.direction = direction;
                s.max_t = direction.norm();
                std::cout << i << "Point " << s.point.x() <<" ," << s.point.y()
                          << "\n Direction:" << s.direction.x() << " ," << s.direction.y() << " Length: " << s.max_t << std::endl;
                sg.add_surface(s);
            }
        }
        sg.write_out_surface();
    }
    // end
    return 0;
}
