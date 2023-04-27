#include "image_converter.h"
/// image converter
/**
 * @fn void imageConverter::read()
 * @brief reads the bmp image into the internal bmp struct
 */
void imageConverter::read() {
    std::filesystem::path file{path};
    if(std::filesystem::exists(file)) {
        std::ifstream bmp_file_input;
        bmp_file_input.open(path,std::ios_base::binary);
        // read file header
        bmp_file_input.read((char*)&bmp.file_header,sizeof(BMPFileHeader_t));
        // check the type
        if(bmp.file_header.file_type != 0x4D42) {
            throw std::runtime_error("Unrecognized bmp file format!");
        }
        // read the info header
        bmp_file_input.read((char*)&bmp.info_header,sizeof(BMPInfoHeader_t));

        // check for additional color headers for transparent images
        if(bmp.info_header.bit_count == 32) {
            // Check if the file has bit mask color information
            if(bmp.info_header.size >= (sizeof(BMPInfoHeader) + sizeof(BMPColorTable_t))) {
                bmp_file_input.read((char*)&bmp.color_header, sizeof(bmp.color_header));
                // Check if the pixel data is stored as BGRA and if the color space type is sRGB
                BMPColorTable_t expected_color_header;
                if(expected_color_header.red_mask != bmp.color_header.red_mask ||
                    expected_color_header.blue_mask != bmp.color_header.blue_mask ||
                    expected_color_header.green_mask != bmp.color_header.green_mask ||
                    expected_color_header.alpha_mask != bmp.color_header.alpha_mask) {
                    throw std::runtime_error("Unexpected color mask format! The program expects the pixel data to be in the BGRA format");
                }
                if(expected_color_header.color_space_type != bmp.color_header.color_space_type) {
                    throw std::runtime_error("Unexpected color space type! The program expects sRGB values");
                }
            } else {
                throw std::runtime_error("Error! Unrecognized file format.");
            }
        }

        // go to the pixel data blocks
        bmp_file_input.seekg(bmp.file_header.offset_data, bmp_file_input.beg);

        // color header is for transparent images
        /// todo prob not necessary
        if(bmp.info_header.bit_count == 32) {
            bmp.info_header.size = sizeof(BMPInfoHeader_t) + sizeof(BMPColorTable_t);
            bmp.file_header.offset_data = sizeof(BMPFileHeader_t) + sizeof(BMPInfoHeader_t) + sizeof(BMPColorTable_t);
        }
        else {
            bmp.info_header.size = sizeof(BMPInfoHeader_t);
            bmp.file_header.offset_data = sizeof(BMPFileHeader_t) + sizeof(BMPInfoHeader_t);
        }

        bmp.file_header.file_size = bmp.file_header.offset_data;


        // resize data for reading
        bmp.data.resize(bmp.info_header.bit_count * bmp.info_header.width * bmp.info_header.height /8);
        // check for row padding
        if(bmp.info_header.width % 4 == 0) {
            // just read the data
            bmp_file_input.read((char*) bmp.data.data(),(long) bmp.data.size());
            bmp.file_header.file_size += static_cast<uint32_t> (bmp.data.size());
        }
        else {
            uint32_t row_stride = bmp.info_header.width * bmp.info_header.bit_count / 8;
            uint32_t new_stride = make_stride_aligned(4,row_stride);
            std::vector<uint8_t> padding_row(new_stride - row_stride);

            for (int y = 0; y < bmp.info_header.height; ++y) {
                bmp_file_input.read((char*)(bmp.data.data() + row_stride * y), row_stride);
                bmp_file_input.read((char*)padding_row.data(), padding_row.size());
            }
            bmp.file_header.file_size += static_cast<uint32_t>(bmp.data.size())
                                         + bmp.info_header.height
                                               * static_cast<uint32_t>(padding_row.size());
        }
    }
    else {
        throw std::invalid_argument("bmp file doesnt exist!!");
    }
}

void imageConverter::detect_colors() {
    // how many datums do we have to combine to get a full input
    int data_format = bmp.info_header.bit_count/8;
    uint32_t full_data = 0;
    int current_shift = 0;
    for(auto part : bmp.data) {
        // add up the data
        full_data |= part << (current_shift * 8);
        // loop controles + saving and comparing
        current_shift++;
        if(current_shift >= data_format) {
            // save data
            compare_save_color_table(full_data);
            // reset
            current_shift = 0;
            full_data = 0;
        }
    }
}

void imageConverter::check_for_white() {
    uint32_t white = WHITE_COLOR_CODE_24_BIT;
    if(!colors_used.contains(white)) {
        std::cerr << "No baseline white for wet nodes detected!" << std::endl;
        /// todo not sure if i should terminate the program here prob not
    }
}

void imageConverter::compare_save_color_table(uint32_t full_color) {
    // easier logic thanks to c++20 :)
    if(!colors_used.contains(full_color)) {
        colors_used.emplace(full_color,full_color);
    }
}

void imageConverter::create_raw() {
    // order the bmp image into a more accessible 2d structure
    // or directly create a 2d struct from the raw bmp data not sure yet
    // create a boundary with the basic sizes can be used in the node generator
    point_t size = {bmp.info_header.width,bmp.info_header.height};
    raw = new rawBoundaryPoints(size);
    point_t current ={0,0};
    // run through raw bmp data
    uint32_t full_data = 0;
    int current_shift = 0;
    int data_format = bmp.info_header.bit_count/8;
    for(auto part : bmp.data) {
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
                // add to the raw boundaries
                raw->read_in_bounce_back(current);
            }
            current = update_position(current);
            // reset
            current_shift = 0;
            full_data = 0;
        }
    }
}

void imageConverter::fill_pkh_with_reformed_raw() {
    pkh.clear();
    for(auto reformed_bp : raw->reformed_boundary_points) {
        pkh.fill_key(reformed_bp->h,reformed_bp->point);
    }
}

void imageConverter::translate_reformed_into_structure() {
    // create a new structure
    current_structure++;
    boundaries->init_structure();
    // generate a pkh association for the reformed nodes
    fill_pkh_with_reformed_raw();
    // init the window to determine weather or not we found a previous partner
    windowedHandles previous_checker(TARGET_SIZE_WINDOW);
    // set up the delete control vector
    std::vector<bool> delete_control(raw->reformed_boundary_points.size(),false);
    // start and current bp
    auto start = raw->reformed_boundary_points.at(0);
    auto current = start;
    // loop over start with first point; x cause we dont need the value
    for(int x = 0; x < raw->reformed_boundary_points.size(); ++x ) {
        // search in all cardinal directions for a partner + add that one
        array_t short_hand = current->point;
        // we use the velocity set for easy directions
        // we could also predict the direction based on the previous one
        /// todo not sure if only cardinal directions should be allowed?!
        for(int i = 1; i < CHANNELS; ++i) {
            // eigen can init a point with an array but cant add stuff up...
            point_t p = short_hand + velocity_set.col(i);
            // search for a viable handle + there always should be one
            handle_t found_handle = pkh.key_translation(p);
            // add logic and move to the next one
            if(found_handle > 0) {
                // check against previous handles only add if a valid handle
                if(!previous_checker.check(found_handle)) {
                    // set up the boundary struct
                    // we intentionally use the wrong handle here
                    boundaries->set_point(current->h,&current->point,current->type);
                    // go to the next element and setup the delete array
                    handle_t array_position = found_handle -1;
                    current = raw->reformed_boundary_points[array_position];
                    delete_control[array_position] = true;
                    // add to the previously found ones
                    previous_checker.add(found_handle);
                    // break out of the for loop
                    break;
                }
            }
            // error
            if(i == CHANNELS-1) {
                throw std::runtime_error("Unclosed surface or long single row!");
            }
        }
        // check weather or not the current point is the start point
        if(current->h == start->h) {
            break;
        }
    }
    /// todo meme analyse prob has more holes than swiss cheese
    // cant delete cause the handle will be invalid after the first run
    // delete/ erase the points in the reformed boundary points now in the structure
    // nessessary to create a seperate flag vector to flag used elements
    // construct a new vector and clear/ delete the old one, move the new one over
    // go over the reformed nodes and add the nodes into raw_boundary points again
    for(int i = 0; i < delete_control.size(); ++i) {
        // only delete not used elements the others get written into the array
        /// important: dont loose pointers to the object in this whole ordeal
        if(delete_control[i]) {
            // delete element
            delete raw->reformed_boundary_points[i];
        }
        else {
            // keep the pointer to the object
            raw->raw_boundary_points.push_back(raw->reformed_boundary_points[i]);
        }
    }
    // clear the reformed array
    raw->reformed_boundary_points.clear();
    // setup the reformed with the the raw data again (move constructor)
    raw->reformed_boundary_points = raw->raw_boundary_points;
}

uint32_t imageConverter::make_stride_aligned(uint32_t align_stride, uint32_t row_stride) {
    uint32_t new_stride = row_stride;
    while (new_stride % align_stride != 0) {
        new_stride++;
    }
    return new_stride;
}

point_t imageConverter::update_position(point_t p) {
    // update the postion based on strides
    point_t size_shorthand = {bmp.info_header.width,bmp.info_header.height};
    p.x()++;
    if(p.x() >= size_shorthand.x()) {
        p.x() = 0;
        p.y()++;
    }
    if(p.y() > size_shorthand.y()) {
        throw std::runtime_error("Image size overrun");
    }
    return p;
}

/// public
imageConverter::imageConverter(std::filesystem::path p) {
    path = p;
}

void imageConverter::init() {
    read();
    detect_colors();
}

void imageConverter::run() {
    raw_run();
    raw_cleanup();
}

void imageConverter::raw_run() {
    // read in the raw data and reduce the set of boundaries considered
    check_for_white();
    create_raw();
    raw->reduce();
}

void imageConverter::raw_cleanup() {
    // delete the temporary raw points work with the reformed points
    raw->delete_raw_boundary_points();
    // start with the reform into boundary structure
    point_t size = {bmp.info_header.width,bmp.info_header.height};
    boundaries = new boundaryPointConstructor(size);
    // init one structure // do while size of reformed boundary is not 0
    translate_reformed_into_structure();
    for(auto bs : boundaries->boundary_structures) {
        bs->rewrite_reformed_boundary_handles();
    }
}

int imageConverter::return_number_of_colors() {
    return colors_used.size();
}

unsigned long imageConverter::return_basic_size() {
    if(bmp.info_header.width > bmp.info_header.height) {
        return  bmp.info_header.width;
    }
    else {
        return bmp.info_header.height;
    }
}