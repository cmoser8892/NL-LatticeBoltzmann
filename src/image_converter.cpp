#include "image_converter.h"

// image converter
/**
 * Reads the bmp image into the internal bmp struct.
 */
void imageConverter::read() {
    bmp_reader->read();
}

/**
 * Goes through the image and lists the colors, not really necessary, more done for statistics.
 */
void imageConverter::detect_colors() {
    // how many datums do we have to combine to get a full input
    int data_format = bmp_reader->bmp.info_header.bit_count/8;
    uint32_t full_data = 0;
    int current_shift = 0;
    for(auto part : bmp_reader->bmp.data) {
        // add up the data
        full_data |= part << (current_shift * 8);
        // loop controles + saving and comparing
        current_shift++;
        if(current_shift >= data_format) {
            // save data
            if(!colors_used.contains(full_data)) {
                colors_used.emplace(full_data,full_data);
            }
            // reset
            current_shift = 0;
            full_data = 0;
        }
    }
}

/**
 * Creates the raw boundary points as given by the image, everything not white (0xfffff) is considered boundary.
 */
void imageConverter::create_raw() {
    // order the bmp image into a more accessible 2d structure
    // or directly create a 2d struct from the raw bmp data not sure yet
    // create a boundary with the basic sizes can be used in the node generator
    point_t size = {bmp_reader->bmp.info_header.width,bmp_reader->bmp.info_header.height};
    raw = new rawBoundaryPoints(size);
    point_t current ={0,0};
    // run through raw bmp data
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
                // add to the raw boundaries
                raw->read_in_bounce_back(current);
            }
            current = update_image_position(current, &size);
            // reset
            current_shift = 0;
            full_data = 0;
        }
    }
}

/**
 * Translates the the reformed raw data into an boundary struct.
 */
void imageConverter::translate_reformed_into_structure() {
    // create a new structure
    current_structure++;
    boundaries->init_structure();
    // generate a pkh association for the reformed nodes
    pkh.clear();
    for(auto reformed_bp : raw->reformed_boundary_points) {
        pkh.fill_key(reformed_bp->h,reformed_bp->point);
    }
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
        // todo not sure if only cardinal directions should be allowed?!
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
    // cant delete cause the handle will be invalid after the first run
    // delete/ erase the points in the reformed boundary points now in the structure
    // nessessary to create a seperate flag vector to flag used elements
    // construct a new vector and clear/ delete the old one, move the new one over
    // go over the reformed nodes and add the nodes into raw_boundary points again
    for(int i = 0; i < delete_control.size(); ++i) {
        // only delete not used elements the others get written into the array
        // important: dont loose pointers to the object in this whole ordeal
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
    raw->rewrite_reformed_boundary_handles();
}

/**
 * Constructor of the image converter.
 * @param p
 */
imageConverter::imageConverter(std::filesystem::path p) {
    bmp_reader = new bmpReader(p);
}

/**
 * Deconstrtor.
 */
imageConverter::~imageConverter() {
    // delete also calls teh deconstructor
    delete raw;
    delete boundaries;
    delete bmp_reader;
}

/**
 * Reads in a bmp image.
 */
void imageConverter::init() {
    read();
    detect_colors();
}

/**
 * Run function will setup everything.
 */
void imageConverter::run() {
    raw_run();
    raw_cleanup();
}

/**
 * Partial run, used for testing.
 */
void imageConverter::raw_run() {
    // read in the raw data and reduce the set of boundaries considered
    if(!check_for_white_wet_nodes()) {
        std::cerr << "No baseline white for wet nodes detected!" << std::endl;
    }
    create_raw();
    raw->reduce();
}

/**
 * Main function to translate.
 */
void imageConverter::raw_cleanup() {
    // delete the temporary raw points work with the reformed points
    raw->delete_raw_boundary_points();
    // start with the reform into boundary structure
    point_t size = {bmp_reader->bmp.info_header.width,bmp_reader->bmp.info_header.height};
    boundaries = new boundaryPointConstructor(size);
    // init one structure // do while size of reformed boundary is not 0
    do {
        // we use the raw boundary points as a helper structure
        raw->raw_boundary_points.clear();
        translate_reformed_into_structure();
    }while(!raw->reformed_boundary_points.empty());
    // boundary handles are the handles from raw data, we have to reform them
    for(auto bs : boundaries->boundary_structures) {
        bs->rewrite_reformed_boundary_handles();
    }
}

/**
 * Return the total number of colors in the bmp image, used for testing.
 * @return the number of colors
 */
int imageConverter::return_number_of_colors() {
    return (int) colors_used.size();
}

/**
 * Check if white is in the image, white is the considered as a wet node.
 * @return yes no white there
 */
bool imageConverter::check_for_white_wet_nodes() {
    return colors_used.contains(WHITE_COLOR_CODE_24_BIT);
}

/**
 * Return either width or height depending which one is bigger.
 * @return returns a long number
 */
unsigned long imageConverter::return_basic_size() {
    if(bmp_reader->bmp.info_header.width > bmp_reader->bmp.info_header.height) {
        return  bmp_reader->bmp.info_header.width;
    }
    else {
        return bmp_reader->bmp.info_header.height;
    }
}