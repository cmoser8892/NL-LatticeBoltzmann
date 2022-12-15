//
// Created by workstation on 16.11.22.
//

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
    NONE = 0,
    DRY,
    WET
}nodeIdentifier_t;

typedef enum boundaryType {
    NO_BOUNDARY = 0, // doesnt do anything equal to an error
    BOUNCE_BACK,
    BOUNCE_BACK_MOVING,
    PRESSURE
}boundaryType_t;
#endif // MY_GENERAL_CODE_TYPES_H
