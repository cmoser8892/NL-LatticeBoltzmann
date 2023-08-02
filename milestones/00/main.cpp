
/*
 * vallgrind call
valgrind --tool=callgrind --dump-instr=yes (p)
 */

#include <opencv2/opencv.hpp>
#include <iostream>

int main() {
    cv::Mat image = cv::imread("/home/sideproject/Downloads/IMG_20221030_131653907_MFNR.jpg");
    if (image.empty()) {
        std::cout << "Error: Unable to load the image." << std::endl;
        return -1;
    }

    cv::imshow("Loaded Image", image);
    cv::waitKey(0);
    return 0;
}