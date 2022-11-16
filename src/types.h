//
// Created by workstation on 16.11.22.
//

#ifndef MY_GENERAL_CODE_TYPES_H
#define MY_GENERAL_CODE_TYPES_H

#include <Eigen/Dense>

using array_t = Eigen::ArrayXd;

typedef enum nodeIdentifier {
    NONE = 0,
    BODY,
    BOUNDARY
}node_identifier_t;

#endif // MY_GENERAL_CODE_TYPES_H
