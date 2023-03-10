#include <iostream>
#include "boundary_point_generator.h"

/**
 * @fn boundaryPointConstructor::boundaryPointConstructor(point_t s)
 * @brief constructor of the boundary points
 * @param s
 */
boundaryPointConstructor::boundaryPointConstructor(point_t s) {
    size = s;
    limits << s.x()-1,s.y()-1;
}

/**
 * @fn boundaryPointConstructor::~boundaryPointConstructor()
 * @brief deconstructor, deletes the elements
 */
boundaryPointConstructor::~boundaryPointConstructor() {
    delete_elements();
}

/**
 * @fn void boundaryPointConstructor::one_direction(int limit, vector_t dir, point_t *start, boundaryType_t  b)
 * @brief constructs a number of boundary points into on direction
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
 * @fn void boundaryPointConstructor::set_point(point_t* p, boundaryType_t b)
 * @brief sets up an individual boundary point
 * @param p
 * @param b
 */
void boundaryPointConstructor::set_point(point_t* p, boundaryType_t b) {
    auto boundary_point = new boundaryPoint_t;
    boundary_point->point = *p;
    boundary_point->type = b;
    boundary_points.push_back(boundary_point);
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
    for(auto b : boundary_points) {
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
    for(auto b : boundary_points) {
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
    for(auto b : boundary_points) {
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
    for(auto b : boundary_points) {
        if(b->point.y() == limit_y) {
            b->type = BOUNCE_BACK_MOVING;
        }
    }
}

/**
 * @fn void boundaryPointConstructor::delete_elements()
 * @brief deletes the elements in the boundary structure
 */
void boundaryPointConstructor::delete_elements() {
    for(auto element : boundary_points) {
        // deletes the element or the pointer
        delete element;
    }
}

void boundaryPointConstructor::visualize_2D_boundary(int size) {
    flowfield_t output;
    output.setZero(size,size);
    for(auto b : boundary_points) {
        output(int(b->point.x()),int(b->point.y())) = 1;
    }
    std::cout << "Boundary-structure" << std::endl;
    std::cout << output << std::endl << std::endl;
}