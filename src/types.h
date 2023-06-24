#ifndef MY_GENERAL_CODE_TYPES_H
#define MY_GENERAL_CODE_TYPES_H

#include <Eigen/Dense>

#define CHANNELS 9 /**< D2Q9 has 9 channels */
#define CARDINAL_DIRECTIONS 4

using array_t = Eigen::ArrayXd; /**< Eigen redefinition, a simple array to hold data */
using matrix_t = Eigen::ArrayXXd; /**< Eigen redefinition, an array with 2 dimensions */
using flowfield_t = Eigen::ArrayXXd; /**< Eigen redefinition, an array with 2 dimensions */

using point_t = Eigen::Vector2d; /**< 2D vector */
using vector_t = Eigen::Vector2d; /**< 2d vector */

using vector3d_t = Eigen::Vector3d; /**< 3d vector */

using handle_t = uint64_t; /**< counter to a node stored in a vector @attention valid handles start with 1 */

using colour_t = uint32_t; /**< placeholder for colour, whatever the right dataformat may be */

extern matrix_t velocity_set; /**< is set as a global */
extern matrix_t weights; /**< is set as a global */
extern matrix_t cardinal_directions; /**< is set as a global */

/**
 * xy (long) coordinates as a struct
 */
typedef struct coordinate {
    long x;
    long y;
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
    PERIODIC,
    PRESSURE_PERIODIC,
    OPEN_INLET,
    OPEN_OUTLET
}boundaryType_t;

/*
 * Note on the toLinks_t:
 * used in the naive implementation and used
 * when constructing neighbourhoods
 */
typedef struct toLinks {
    int channel;
    handle_t handle; // form 1 to n -1 for valid handles
}toLinks_t;

typedef struct simulation_parameters {
    double relaxation = 0.5;
    double u_wall = 0;
    double dt = 1;
}simulation_parameters_t;

// eigen smart pointer
using link_pointer = Eigen::internal::pointer_based_stl_iterator<array_t>;
using array_pointer = link_pointer;
#endif // MY_GENERAL_CODE_TYPES_H
