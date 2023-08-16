
#include <opencv2/opencv.hpp>
#include "helper_functions.h"
int main() {
    auto test_image = get_base_path();
    test_image.append("tests");
    test_image.append("test_images");
    test_image.append("wiggle_lines.png");
    // call the drawer
    // Load an image using OpenCV
    cv::Mat image = cv::imread(test_image.string()); // Replace with the path to your image

    // Convert the image to grayscale
    cv::Mat grayImage;
    cv::cvtColor(image, grayImage, cv::COLOR_BGR2GRAY);

    // Apply Canny edge detection to the grayscale image with higher thresholds
    cv::Mat edges;
    cv::Canny(grayImage, edges, 100, 200); // Increased thresholds

    // Find contours in the edges image
    std::vector<std::vector<cv::Point>> contours;
    cv::findContours(edges, contours, cv::RETR_EXTERNAL, cv::CHAIN_APPROX_SIMPLE);

    // Fit lines to the detected contours using least squares with higher minimum points
    std::vector<cv::Vec4f> lines;
    for (const auto& contour : contours) {
        if (contour.size() > 20) { // Increased minimum points
            cv::Vec4f line;
            cv::fitLine(contour, line, cv::DIST_L2, 0, 0.01, 0.01);

            // Draw multiple short lines along the detected line
            int segmentLength = 20; // Length of each line segment
            for (int i = 0; i < contour.size() - segmentLength; i += segmentLength) {
                cv::Point pt1(line[2] + i * line[0], line[3] + i * line[1]);
                cv::Point pt2(line[2] + (i + segmentLength) * line[0], line[3] + (i + segmentLength) * line[1]);
                cv::line(image, pt1, pt2, cv::Scalar(0, 0, 255), 2, cv::LINE_AA);
            }
        }
    }

    // Display the image with detected lines
    cv::imshow("Image with Detected Lines", image);
    cv::waitKey(0);
}