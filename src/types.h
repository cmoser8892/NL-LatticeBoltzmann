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
 * In the final implementation it should be handle + channel to directly
 * would need extra processing to get the combo thou
 * access the correct node and channel.
 * To get the channel: handle % Number_of_channels
 * To get the correct node: handle / Number_of_channels
 * maybe some shenanigans with 0 and max size_t
 * currently handles are only valid if > 0 and smaller < Max size
 * may have to be changed, also we get like 9e18 valid nodes
 */
typedef struct toLinks {
    int channel;
    handle_t handle; // form 1 to n -1 for valid handles
}toLinks_t;

#endif // MY_GENERAL_CODE_TYPES_H
