#include "simulation.h"
#include "image_converter.h"
#include <iostream>
#include <chrono>
/**
 * 10 main, imb integration, we simulate a quader of a fluid rotating, boarders are off.
 * @param argc
 * @param argv
 * @return
 */
int main(int argc, char *argv[]) {
    point_t starter = {3.1,3.1};
    straight_t input;
    straightGenerator sg;
    long canvas_size = 100;
    double side_length = 69; // with a distance of 0.75 we should get 80 markers
    // we put in a quader
    input.point = starter;
    input.direction = {0,1};
    input.max_t = side_length;
    sg.add_surface(input);
    input.point += input.direction * side_length;
    input.direction = {1,0};
    input.max_t = side_length;
    sg.add_surface(input);
    input.point += input.direction * side_length;
    input.direction = {0,-1};
    input.max_t = side_length;
    sg.add_surface(input);
    input.point += input.direction * side_length;
    input.direction = {-1,0};
    input.max_t = side_length;
    sg.add_surface(input);
    sg.surface_mass_center();
    nodeGenerator ng(&sg);
    ng.init_surface(canvas_size);
    ng.visualize_2D_nodes();
    std::cout << ng.node_infos.size();
    return 0;
}