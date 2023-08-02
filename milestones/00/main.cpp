
/*
 * vallgrind call
valgrind --tool=callgrind --dump-instr=yes (p)
 */

#include <iostream>
#include <vector>
#include <opencv2/opencv.hpp>
#include "helper_functions.h"
using namespace cv;
using namespace std;

int main() {
    // Load the image
    auto test_image = get_base_path();
    test_image.append("tests");
    test_image.append("test_images");
    test_image.append("cool_duck.png");
    Mat image = imread(test_image.string(), IMREAD_GRAYSCALE);
    if (image.empty()) {
        cout << "Error: Unable to load the image." << endl;
        return -1;
    }

    // Parameters for Harris Corner Detection
    int blockSize = 2; // Size of the neighborhood considered for corner detection
    int apertureSize = 3; // Aperture parameter for Sobel operator
    double k = 0.04; // Harris corner free parameter
    double threshold = 0.01; // Threshold to distinguish corners from non-corners

    // Detect corners using Harris Corner Detection
    Mat cornerStrength;
    cornerHarris(image, cornerStrength, blockSize, apertureSize, k);

    // Compute the maximum corner strength value
    double maxStrength;
    minMaxLoc(cornerStrength, nullptr, &maxStrength);

    // Create a binary image to store the corner locations
    Mat cornerBinary = Mat::zeros(cornerStrength.size(), CV_8U);

    // Iterate through the corner strength matrix and set the binary image pixels
    // to 255 if the corresponding corner strength is above the threshold
    for (int y = 0; y < cornerStrength.rows; ++y) {
        for (int x = 0; x < cornerStrength.cols; ++x) {
            if (cornerStrength.at<float>(y, x) > threshold * maxStrength) {
                cornerBinary.at<uchar>(y, x) = 255;
            }
        }
    }

    // Find the non-zero (corner) locations
    vector<Point> cornerPoints;
    cout << cornerPoints.size() << std::endl;
    findNonZero(cornerBinary, cornerPoints);
    // Display the detected corners
    if(0) {
        // Draw circles around the detected corners
        for (const Point& p : cornerPoints) {
            circle(image, p, 5, Scalar(255), 2);
        }
        imshow("Detected Corners", image);
        waitKey(0);
    }

    // Compute the convex hull of the detected corners
    vector<Point> convexHullPoints;
    convexHull(cornerPoints, convexHullPoints);

    // Draw the convex hull on the original image
    vector<vector<Point>> hulls = {convexHullPoints};
    drawContours(image, hulls, 0, Scalar(255), 2);

    // Display the image with the detected corners and the recovered surface
    imshow("Detected Corners and Surface", image);
    string output = "out.jpg";
    imwrite(output,image);
    waitKey(0);


    return 0;
}
