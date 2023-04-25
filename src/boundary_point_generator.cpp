#include <iostream>
#include "boundary_point_generator.h"

rawBoundaryPoints::rawBoundaryPoints(point_t s) {
    size = s;
}

rawBoundaryPoints::~rawBoundaryPoints() {
    delete_raw_boundary_points();
    delete_reformed_boundary_points();
}

/**
 * @fn void rawBoundaryPoints::delete_raw_boundary_points()
 * @brief deletes the boundary point structs saved in the raw boundary vector
 * @attention clearing the vector prevents double deletion of dangeling pointers
 */
void rawBoundaryPoints::delete_raw_boundary_points() {
    for(auto d : raw_boundary_points) {
        delete d;
    }
    raw_boundary_points.clear();
}

/**
 * @fn void rawBoundaryPoints::delete_reformed_boundary_points()
 * @brief deletes the boundary point structs saved in the reformed boundary vector
 * @attention clearing the vector prevents double deletion of dangeling pointers
 */
void rawBoundaryPoints::delete_reformed_boundary_points() {
    for(auto d : reformed_boundary_points) {
        delete d;
    }
    reformed_boundary_points.clear();
}

/**
 * @fn void rawBoundaryPoints::fill_keys()
 * @brief for loop over rbps to fill the pointkeyhash class
 */
void rawBoundaryPoints::fill_keys() {
    for(auto rbp : raw_boundary_points) {
        pkh.fill_key(rbp->h,rbp->point);
    }
}


void rawBoundaryPoints::visualize_2D_boundary() {
    flowfield_t output;
    output.setZero(size.x(),size.y());
    for(auto b : raw_boundary_points) {
        ++output(int(b->point.x()),int(b->point.y()));
    }
    std::cout << "Boundary-structure" << std::endl;
    std::cout << output << std::endl << std::endl;
}

void rawBoundaryPoints::rewrite_reformed_boundary_handles() {
    handle_t start = 0;
    for(auto reformed_bp : reformed_boundary_points) {
        reformed_bp->h = ++start;
    }
}

rawBoundaryPoints::border_return_code_t rawBoundaryPoints::check_boarder(boundaryPoint_t &b) {
    border_return_code_t return_code = INSIDE;
    point_t short_hand = b.point;
    int counter = 0;
    //
    if(short_hand.x() == 0) {
        counter++;
    }
    if(short_hand.y() == 0) {
        counter++;
    }
    if(short_hand.x() == size.x()) {
        counter++;
    }
    if(short_hand.y() == size.y()) {
        counter++;
    }
    // assigning return codes
    if(counter > 2) {
        throw std::runtime_error("Undefined point, not possible in 2D");
    }
    else if(counter == 2) {
        return_code = CORNER;
    }
    else if(counter == 1) {
        return_code = BOARDER;
    }
    else {
        // nop
    }
    return return_code;
}

int rawBoundaryPoints::set_max_neighbors(rawBoundaryPoints::border_return_code_t b) {
    int return_code = -1;
    switch (b) {
    case CORNER:
        return_code = 3;
        break;
    case BOARDER:
        return_code = 5;
        break;
    case INSIDE:
        return_code = 8;
        break;
    default:
        throw std::runtime_error("Default in switch case");
    }
    return return_code;
}

void rawBoundaryPoints::read_in_bounce_back(point_t p) {
    coordinate_t coordinate;
    coordinate.x = std::floor(p.x());
    coordinate.y = std::floor(p.y());
    read_in_bounce_back(coordinate);
}

void rawBoundaryPoints::read_in_bounce_back(coordinate_t coordinate) {
    auto new_bp = new boundaryPoint_t;
    new_bp->h = ++current_handle;
    // aka should have done it the other way around
    new_bp->point.x() = (double) coordinate.x;
    new_bp->point.y() = (double) coordinate.y;
    // setup identification variables
    new_bp->dw = DRY;
    new_bp->type = BOUNCE_BACK;
    // push into structure
    raw_boundary_points.push_back(new_bp);
}

void rawBoundaryPoints::reduce() {
    /// todo exceptionally similar to determine neighbors in Neighborhood
    /// todo still quite a number of linear searches
    // gets rid of all the unnecessary boundary points
    // build hash search space
    fill_keys();
    handle_t start = 0;
    // loop over the raw nodes
    for(auto rbp : raw_boundary_points) {
        // determine boarder code
        border_return_code_t cases = check_boarder(*rbp);
        // set the number on set case where we discard a point
        int out_number = set_max_neighbors(cases);
        // go in all directions to look for a neighbor
        int found_neighbors = 0;
        array_t short_hand = rbp->point;
        // we use the velocity set for easy directions
        for(int i = 1; i < CHANNELS; ++i) {
            // eigen can init a point with an array but cant add stuff up...
            point_t current = short_hand + velocity_set.col(i);
            handle_t found_handle = pkh.key_translation(current);
            if(found_handle > 0) {
                found_neighbors++;
            }

        }
        // check how many neighbors were found
        if(found_neighbors >= out_number) {
            // out nop
        }
        else {
            // add to the reformed nodes
            auto copy  = new boundaryPoint_t;
            // copy over
            copy->h = ++start;
            copy->point = rbp->point;
            copy->dw = rbp->dw;
            copy->type = rbp->type;
            // add to the vector
            reformed_boundary_points.push_back(copy);
        }
    }
    // clear the hash table for later use
    pkh.clear();
}


/**
 * @fn boundaryStructure::~boundaryStructure()
 * @brief de-construcutur deletes all the boundary points
 */
boundaryStructure::~boundaryStructure() {
    for(auto d : boundary_points) {
        delete d;
    }
    // not really necessary here but good practise
    boundary_points.clear();
}

/**
 * @fn boundaryPointConstructor::boundaryPointConstructor(point_t s)
 * @brief constructor of the boundary points
 * @param s
 */
boundaryPointConstructor::boundaryPointConstructor(point_t s) {
    size = s;
    limits << s.x()-1,s.y()-1;
    current_structure = -1;

}

/**
 * @fn boundaryPointConstructor::~boundaryPointConstructor()
 * @brief deconstructor, deletes the elements
 */
boundaryPointConstructor::~boundaryPointConstructor() {
    delete_structures();
}
/**
 * @fn void boundaryPointConstructor::init_structure()
 * @brief initilizes a boundary structure to better describe disjunct boundariess
 */
void boundaryPointConstructor::init_structure() {
    auto bs = new boundaryStructure;
    added_handle = 0;
    boundary_structures.push_back(bs);
    current_structure++;

}

/**
 * @fn void boundaryPointConstructor::one_direction(int limit, vector_t dir, point_t *start, boundaryType_t  b)
 * @brief constructs a number of boundary points into on direction only really good for cardinal directions
 * @param limit
 * @param dir
 * @param start
 * @param b
 */
void boundaryPointConstructor::one_direction(int limit, vector_t dir, point_t *start, boundaryType_t  b) {
    // goes up and down a line and inizlizes it
    for(int i = 0; i < limit; i++) {
        set_point(start,b);
        *start += dir;
    }
}

/**
 * @fn void boundaryPointConstructor::steps_direction(int steps, vector_t dir, point_t *start, boundaryType_t b)
 * @brief creates step in one direction (similar to one direction
 * @param steps
 * @param dir
 * @param start
 * @param b
 */
void boundaryPointConstructor::steps_direction(int steps, vector_t dir, point_t *start, boundaryType_t b) {
    // creation of x corners
    for(int i = 0; i < steps; ++i) {
        corner_creation(dir,start,b);
    }
}

/**
 * @fn void boundaryPointConstructor::corner_creation(vector_t dir, point_t *start, boundaryType_t b)
 * @brief creates a corner
 */
void boundaryPointConstructor::corner_creation(vector_t dir, point_t *start, boundaryType_t b) {
    // creation of a corner
    vector_t normal = {dir.y(), -dir.x()};
    set_point(start,b);
    vector_t addup = 0.5 *(normal+dir);
    *start += addup;
    set_point(start,b);
    addup = 0.5 *(-normal + dir);
    *start += addup;
}
/**
 * @fn void boundaryPointConstructor::set_point(point_t* p, boundaryType_t b)
 * @brief sets up an individual boundary point
 * @param p
 * @param b
 */
void boundaryPointConstructor::set_point(point_t* p, boundaryType_t b) {
    set_point(++added_handle,p,b);
}


void boundaryPointConstructor::set_point(handle_t h, point_t *p, boundaryType_t b) {
    if(current_structure<0) {
        std::cerr << "No structure" << std::endl;
        return;
    }
    // gerenate the boundary node
    auto boundary_point = new boundaryPoint_t;
    boundary_point->h = h;
    boundary_point->point = *p;
    boundary_point->dw = DRY;
    boundary_point->type = b;
    boundary_structures.at(current_structure)->boundary_points.push_back(boundary_point);
}

void boundaryPointConstructor::rewrite_handles() {
    for(auto bs : boundary_structures) {
        handle_t start = 0;
        for (auto bp : bs->boundary_points) {
            bp->h = ++start;
        }
    }
}
/**
 * @fn void boundaryPointConstructor::init_quader()
 * @brief sets up a quader of boundary points
 */
void boundaryPointConstructor::init_quader() {
    point_t current;
    current.setZero();
    init_quader(current);
}

/**
 * @fn void boundaryPointConstructor::init_chopped_quader(point_t point, int devider)
 * @brief sets up a quader where are are part was chopped off
 * @param point
 * @param devider
 */
void boundaryPointConstructor::init_chopped_quader(point_t point, int devider) {
    // devider is a devider
    if(devider == 0) {
        devider = 2147483647;
    }
    if(devider < 2) {
        throw std::runtime_error("unrealitic devider");
    }
    init_structure();
    boundaryType_t type = BOUNCE_BACK;
    point_t current = point;
    // go from 0 till the in x directions
    int chop_x = int(size.x()/ devider);
    int chop_y = int(size.y()/ devider);
    vector_t direction;
    // go through x
    direction = {1,0};
    one_direction(int(limits.x())-chop_x,direction,&current, type);
    direction = {0,1};
    one_direction(chop_y,direction,&current, type);
    direction = {1,0};
    one_direction(chop_x,direction,&current, type);
    direction = {0,1};
    one_direction(int(limits.y())-chop_y,direction,&current, type);
    // go through y
    direction = {-1,0};
    one_direction(int(limits.x()),direction,&current, type);
    // go through x
    direction = {0,-1};
    one_direction(int(limits.y()),direction,&current, type);
}

/**
 * @fn void boundaryPointConstructor::init_quader(point_t p)
 * @brief sets up a quader that doesnt start at 0,0
 * @param p
 */
void boundaryPointConstructor::init_quader(point_t p) {
    init_quader(p,size);
}

/**
 * @fn void boundaryPointConstructor::init_quader(point_t p,vector_t size)
 * @brief sets up a quader with a specific size
 * @param p
 * @param size
 */
void boundaryPointConstructor::init_quader(point_t p,vector_t size) {
    boundaryType_t type = BOUNCE_BACK;
    point_t current = p;
    init_structure();
    // go from 0 till the in x directions
    vector_t direction;
    int size_x = int(size.x()-1);
    int size_y = int(size.y()-1);
    // go through x
    direction = {1,0};
    one_direction(size_x,direction,&current, type);
    // go through x
    direction = {0,1};
    one_direction(size_y,direction,&current, type);
    // go through x
    direction = {-1,0};
    one_direction(size_x,direction,&current, type);
    // go through x
    direction = {0,-1};
    one_direction(size_y,direction,&current, type);
}

/**
 * @fn void boundaryPointConstructor::init_sliding_lid()
 * @brief sets up a sliding lid
 */
void boundaryPointConstructor::init_sliding_lid() {
    // greate a slinding lid container with the given sizes
    //in our case y max is the boundary that is moving
    // so all boundaries with y = y_max are BOUNDARY_MOVING
    // init a quader
    double limit_y = limits.y();
    init_quader();
    for(auto b : boundary_structures.at(0)->boundary_points) {
        if(b->point.y() == limit_y) {
            b->type = BOUNCE_BACK_MOVING;
        }
    }
}

/**
 * @fn void boundaryPointConstructor::init_chopped_sliding_lid(point_t start, int chopfactor)
 * @brief sets up a sliding lid where a part is missing
 * @param start
 * @param chopfactor
 */
void boundaryPointConstructor::init_chopped_sliding_lid(point_t start, int chopfactor) {
    double limit_y = limits.y() + start.y();
    init_chopped_quader(start,chopfactor);
    for(auto b : boundary_structures.at(0)->boundary_points) {
        if(b->point.y() == limit_y) {
            b->type = BOUNCE_BACK_MOVING;
        }
    }
}

/**
 * @fn void boundaryPointConstructor::init_quader_side_chopped(point_t start, int chopsize)
 * @brief quader with part missing can be anywhere
 * @param start
 * @param chopsize
 */
void boundaryPointConstructor::init_quader_side_chopped(point_t start, int chopsize) {
    if(chopsize >= size.x()) {
        throw std::runtime_error("recheck sizes");
    }
    if(chopsize >= size.y()) {
        throw std::runtime_error("recheck sizes");
    }
    boundaryType_t type = BOUNCE_BACK;
    point_t current = start;
    int size_x = int(limits.x());
    int size_y = int(limits.y());
    // go from 0 till the in x directions
    vector_t direction;
    // go through x
    direction = {1,0};
    one_direction(size_x,direction,&current, type);
    // go through y side
    direction = {0,1};
    one_direction((size_y-chopsize)/2,direction,&current, type);
    // cut
    direction = {-1,0};
    one_direction(chopsize,direction,&current,type);
    direction = {0,1};
    one_direction(chopsize,direction,&current,type);
    direction = {1,0};
    one_direction(chopsize,direction,&current,type);
    // fill the rest of the side
    direction = {0,1};
    one_direction((size_y-chopsize)/2 + (size_y-chopsize)%2,direction,&current, type);
    // go through y
    direction = {-1,0};
    one_direction(int(limits.x()),direction,&current, type);
    // go through x
    direction = {0,-1};
    one_direction(int(limits.y()),direction,&current, type);
}

/**
 * @fn void boundaryPointConstructor::init_sliding_lid_side_chopped(point_t start, int chopsize)
 * @brief choped of sliding lid
 * @param start
 * @param chopsize
 */
void boundaryPointConstructor::init_sliding_lid_side_chopped(point_t start, int chopsize) {
    double limit_y = limits.y() + start.y();
    init_quader_side_chopped(start,chopsize);
    for(auto b : boundary_structures.at(0)->boundary_points) {
        if(b->point.y() == limit_y) {
            b->type = BOUNCE_BACK_MOVING;
        }
    }
}

/**
 * @fn void boundaryPointConstructor::init_sliding_lid_inner(point_t start, point_t continues, vector_t inner_size)
 * @brief sliding lid with an quadratic object in it
 * @param start
 * @param continues
 * @param inner_size
 */
void boundaryPointConstructor::init_sliding_lid_inner(point_t start, point_t continues, vector_t inner_size) {
    init_quader(start,size);
    init_quader(continues,inner_size);
    double limit_y = limits.y() + start.y();
    for(auto b : boundary_structures.at(0)->boundary_points) {
        if(b->point.y() == limit_y) {
            b->type = BOUNCE_BACK_MOVING;
        }
    }
}

/**
 * @fn void boundaryPointConstructor::delete_structures()
 * @brief deletes the structures
 */
void boundaryPointConstructor::delete_structures() {
    for(auto bs: boundary_structures) {
        delete bs;
    }
    boundary_structures.clear();
    current_structure = -1;
}

/**
 * @fn void boundaryPointConstructor::visualize_2D_boundary(int size)
 * @brief simple visualizer for boundaries
 * @param size
 */
void boundaryPointConstructor::visualize_2D_boundary(int size) {
    flowfield_t output;
    output.setZero(size,size);
    for(auto bs : boundary_structures) {
        for(auto b : bs->boundary_points) {
            ++output(int(b->point.x()),int(b->point.y()));
        }
    }
    std::cout << "Boundary-structure" << std::endl;
    std::cout << output << std::endl << std::endl;
}

/**
 * @fn long boundaryPointConstructor::total_boundary_nodes()
 * @brief sums up all the boundary points
 * @return sum of all the boundary points, indipendent of the individual stuctures formed
 */
long boundaryPointConstructor::total_boundary_nodes() {
    long number_of_nodes = 0;
    for(auto const bs : boundary_structures) {
        number_of_nodes += long(bs->boundary_points.size());
    }
    return number_of_nodes;
}

/**
 * @fn void boundaryPointConstructor::init_poiseuille_flow()
 * @brief standard p flow constructor
 */
void boundaryPointConstructor::init_poiseuille_flow() {
    // init a quader and relable the sides
    init_quader();
    // 0 side
    double limit_x = 0;
    for(auto b : boundary_structures.at(0)->boundary_points) {
        if(b->point.x() == limit_x) {
            b->dw = WET;
            b->type = PRESSURE_PERIODIC;
        }
    }
    // maxsize side
    limit_x = limits.x();
    for(auto b : boundary_structures.at(0)->boundary_points) {
        if(b->point.x() == limit_x) {
            b->dw = WET;
            b->type = PRESSURE_PERIODIC;
        }
    }
    // relabel the top and bottom layer
    double limit_y = 0;
    for(auto b : boundary_structures.at(0)->boundary_points) {
        if(b->point.y() == limit_y) {
            b->dw = DRY;
            b->type = BOUNCE_BACK;
        }
    }
    limit_y = limits.y();
    for(auto b : boundary_structures.at(0)->boundary_points) {
        if(b->point.y() == limit_y) {
            b->dw = DRY;
            b->type = BOUNCE_BACK;
        }
    }
    // nessessary to sort him
    std::sort(boundary_structures[0]->boundary_points.begin(), boundary_structures[0]->boundary_points.end(),  sorter_wet_dry_boundaries);
}

void boundaryPointConstructor::pressure_inlet() {
    // init a quader and relable the sides
    init_quader();
    // 0 side
    double limit_x = 0;
    for(auto b : boundary_structures.at(0)->boundary_points) {
        if(b->point.x() == limit_x) {
            b->dw = DRY;
            b->type = OPEN_INLET;
        }
    }
    // maxsize side
    limit_x = limits.x();
    for(auto b : boundary_structures.at(0)->boundary_points) {
        if(b->point.x() == limit_x) {
            b->dw = DRY;
            b->type = BOUNCE_BACK;
        }
    }
    // relabel the top and bottom layer
    double limit_y = 0;
    for(auto b : boundary_structures.at(0)->boundary_points) {
        if(b->point.y() == limit_y) {
            b->dw = DRY;
            b->type = BOUNCE_BACK;
        }
    }
    limit_y = limits.y();
    for(auto b : boundary_structures.at(0)->boundary_points) {
        if(b->point.y() == limit_y) {
            b->dw = DRY;
            b->type = BOUNCE_BACK;
        }
    }
    // nessessary to sort him
    std::sort(boundary_structures[0]->boundary_points.begin(), boundary_structures[0]->boundary_points.end(),  sorter_wet_dry_boundaries);
}

/**
 * @fn bool sorter_wet_dry_boundaries(boundaryPoint_t * p1, boundaryPoint_t * p2)
 * @brief function that compares the wet and dry state of a node, wet -> 2, dry -> 1
 * @param p1
 * @param p2
 * @return true/false
 */
bool sorter_wet_dry_boundaries(boundaryPoint_t * p1, boundaryPoint_t * p2) {
    // Wet = 2, Dry = 1
    return (p1->dw > p2->dw);
}