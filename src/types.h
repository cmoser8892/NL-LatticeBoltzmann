//
// Created by workstation on 16.11.22.
//

#ifndef MY_GENERAL_CODE_TYPES_H
#define MY_GENERAL_CODE_TYPES_H

#include <Eigen/Dense>

using array_t = Eigen::ArrayXd;
using matrix_t = Eigen::ArrayXXd;

array_t velocity_d2q9 =  {{1,1,0,-1,0 ,1,-1,-1, 1}
                        ,{0,0,1,0 ,-1,1, 1,-1,-1}};

typedef enum nodeIdentifier {
    NONE = 0,
    BODY,
    BOUNDARY
}node_identifier_t;

#endif // MY_GENERAL_CODE_TYPES_H
