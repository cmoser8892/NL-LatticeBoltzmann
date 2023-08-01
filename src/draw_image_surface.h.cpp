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
}

void surfaceDrawer::run() {

}