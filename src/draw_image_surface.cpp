#include "drawn_image_surface.h"

// private
void surfaceDrawer::read() {
    bmp_reader->read();
}

void surfaceDrawer::convert_points() {
    // read in the points
    point_t size = {bmp_reader->bmp.info_header.width,bmp_reader->bmp.info_header.height};
    point_t current ={0,0};
    // read into data structure
    uint32_t full_data = 0;
    int current_shift = 0;
    int data_format = bmp_reader->bmp.info_header.bit_count/8;
    for(auto part : bmp_reader->bmp.data) {
        // add up the data
        full_data |= part << (current_shift * 8);
        // loop controles + saving and comparing
        current_shift++;
        if(current_shift >= data_format) {
            // save data and update the position
            if(full_data == WHITE_COLOR_CODE_24_BIT) {
                // nop
            }
            else {
                // save the point somewhere
                points.push_back(current);
            }
            current = update_image_position(current, &size);
            // reset
            current_shift = 0;
            full_data = 0;
        }
    }
}

void surfaceDrawer::fill_hashtable() {
    for(auto point : points) {
        rpkh.fill_key(handle_runner,point);
        handle_runner++;
    }
 }

 point_t surfaceDrawer::interpolate_around(point_t p) {
    std::vector<handle_t> relevant_points = rpkh.ranging_key_translation(p,range);
    point_t mid = p;
    for(auto h : relevant_points) {
        mid += points[h];
    }
    mid /= ((double)relevant_points.size() +1);
    return mid;
 }

 vector_t surfaceDrawer::determine_init_surface_direction(point_t p, double jump) {
    // go in all major directions and note the one where we find the most results
    std::vector<unsigned long> sizes;
    for(int k = 0; k < major_directions.cols(); ++k) {
        vector_t current = p + point_t(major_directions.col(k))*jump;
        std::vector<handle_t> candidates = rpkh.ranging_key_translation(current,range);
        sizes.push_back(candidates.size());
    }
    // look for the on were we found the most
    int maxElementIndex = (int)std::distance(sizes.begin(),std::max_element(sizes.begin(),sizes.end()));
    vector_t initial = major_directions.col(maxElementIndex);
    return initial;
 }

 bool surfaceDrawer::look_for_last(point_t current) {
    handle_t knock_it_off = points.size()-1;
    return rpkh.ranging_key_look_for_specific(current,range,knock_it_off);
 }

 void surfaceDrawer::add_surface(point_t current, point_t previous) {
     vector_t direction = current - previous;
     straight_t  s;
     s.point = previous;
     s.direction = direction;
     s.max_t = direction.norm();
     surface_storage.add_surface(s);
 }

// public
surfaceDrawer::surfaceDrawer(std::filesystem::path p) {
    bmp_reader = new bmpReader(p);
}

surfaceDrawer::~surfaceDrawer() {
    delete bmp_reader;
}

void surfaceDrawer::init() {
    read();
    convert_points();
    fill_hashtable();
}

void surfaceDrawer::run() {
    point_t current = points[0];
    int total_steps = int((double)points.size() * 2 * 1/step);
    // determine initials
    current = interpolate_around(current);
    point_t previous = current;
    // add the first point last to the points vector
    points.push_back(current);
    rpkh.fill_key(handle_runner,current);
    vector_t surface_direction = determine_init_surface_direction(current,1);
    // for loop
    for(int i = 0; i < total_steps; ++i) {
        // go step in the old direction
        current += surface_direction.normalized()*step;
        // look around and correct
        current = interpolate_around(current);
        // add as a surface
        add_surface(current,previous);
        std::cout << current.x() << " " << current.y() << std::endl;
        std::cout << surface_direction.x() << " " << surface_direction.y() << std::endl;
        // switch over
        surface_direction  = (current - previous).normalized();
        previous = current;
        // every so often we look further
        // surface_direction = determine_init_surface_direction(current,3);
        // check for the end point
        if(look_for_last(current)) {
            break;
        }
    }
    // draw the last surface
    add_surface(points[points.size()-1],previous);
}