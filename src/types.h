#ifndef MY_GENERAL_CODE_TYPES_H
#define MY_GENERAL_CODE_TYPES_H

#include <Eigen/Dense>

#define CHANNELS 9

// redefinitions
using array_t = Eigen::ArrayXd;
using matrix_t = Eigen::ArrayXXd;
using flowfield_t = Eigen::ArrayXXd;

using point_t = Eigen::Vector2d;
using vector_t = Eigen::Vector2d;

using handle_t = uint64_t; // counter to denote a node

// globals (defined in functions)
extern matrix_t velocity_set;
extern matrix_t weights;

// structs
typedef struct coordinate {
    unsigned x;
    unsigned y;
}coordinate_t;

typedef enum nodeIdentifier {
    UNKNOWN = 0,
    DRY,
    WET
}nodeIdentifier_t;

typedef enum boundaryType {
    NO_BOUNDARY = 0, // doesnt do anything equal to an error
    BOUNCE_BACK,
    BOUNCE_BACK_MOVING,
    PRESSURE_PERIODIC
}boundaryType_t;

/**
 * Note on the toLinks_t:
 * used in the naive implementation and used
 * when constructing neighbourhoods
 */
typedef struct toLinks {
    int channel;
    handle_t handle; // form 1 to n -1 for valid handles
}toLinks_t;

// eigen smart pointer
using link_pointer = Eigen::internal::pointer_based_stl_iterator<array_t>;

#endif // MY_GENERAL_CODE_TYPES_H
