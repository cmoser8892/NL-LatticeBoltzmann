/**
 * Thoughts:
 * got more or less two possiblities to implement moving bbs
 * do it addhoc style
 * or create a full boundary representation of a sliding lid and then use that
 */
#ifndef MY_GENERAL_CODE_TYPES_H
#define MY_GENERAL_CODE_TYPES_H

#include <Eigen/Dense>

#define CHANNELS 9

using array_t = Eigen::ArrayXd;
using matrix_t = Eigen::ArrayXXd;
using flowfield_t = Eigen::ArrayXXd;

using point_t = Eigen::Vector2d;
using vector_t = Eigen::Vector2d;

using handle_t = size_t; // not 0 and max_size
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
 * In the final implemntation it should be handle + channel to directlly
 * access the correct node and channel.
 * To get the channel: handle % Number_of_channels
 * To get the correct node: handle / Number_of_channels
 * maybe some shninanigans with 0 and max size_t
 * currentlly handles are only valid if > 0 and smaller < Max size
 * may have to be changed, also we get like 9e18 valid nodes
 */
typedef struct toLinks {
    int channel;
    handle_t handle; // form 1 to n -1 for valid handles
}toLinks_t;

#endif // MY_GENERAL_CODE_TYPES_H
