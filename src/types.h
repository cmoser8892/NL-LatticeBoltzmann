#ifndef MY_GENERAL_CODE_TYPES_H
#define MY_GENERAL_CODE_TYPES_H

#include <Eigen/Dense>

#define CHANNELS 9 /**< D2Q9 has 9 channels */
#define CARDINAL_DIRECTIONS 4
#define KERNEL_3_FACTOR 2 /**< When we give a range for IBM_OUTER we have to multiply with 2 to get the correct cutoff */

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
extern matrix_t major_directions; /**< is set as a global */

/**
 * xy (long) coordinates as a struct
 */
typedef struct coordinate {
    long x; /**<  well x */
    long y; /**<  y well */
}coordinate_t;

/**
 * Dry Wet enum identifier of nodes
 */
typedef enum nodeIdentifier {
    UNKNOWN = 0, /**< Neither wet nor dry (undifined) */
    DRY, ///
    WET
}nodeIdentifier_t;

/**
 * Type of Boundary
 * @note the existence of a tag does not guaranty that there is actual an implementation
 */
typedef enum boundaryType {
    INIT_NONE = 0,
    NO_BOUNDARY, /**<  Just a fluid node */
    BOUNCE_BACK, /**< Bounce back type boundary */
    BOUNCE_BACK_MOVING, /**<  Moving bounce back type boundary */
    PERIODIC, /**<  periodic type boundary */
    PRESSURE_PERIODIC, /**<  pressure periodic tag */
    IBM_OUTER,          /**< Immersed Boundary Method tag, outside of the boundary markers */
    IBM_INNER,          /**< Immersed Boundary Method tag, inside of the boundary markers */
    FAKE_FORCEING_INNER,
    FAKE_FORCEING,
    OPEN_INLET,
    OPEN_OUTLET
}boundaryType_t;

/**
 * Kernel identifiers used for LBM boundaries.
 * @ref Viggen P.469
 */
typedef enum kernelType {
    KERNEL_A = 0, /**< Phi_2(x) -> we need look 2 lattice sites around the marker */
    KERNEL_B,     /**< Phi_3(x) -> we need look 3 lattice sites around the marker */
    KERNEL_C      /**< Phi_4(x) -> we need look 4 lattice sites around the marker */
}kernelType_t;

/**
 * Link to other neighboring nodes
 * @note on the toLinks_t: used in the naive implementation and used when constructing neighbourhoods
 */
typedef struct toLinks {
    int channel; /**<  channel to be streamed to */
    handle_t handle; /**<  form 1 to n -1 for valid handles */
}toLinks_t;

/**
 * Simulation parameter struct.
 */
typedef struct simulation_parameters {
    double relaxation = 0.5; /**<  current relaxation */
    double u_wall = 0; /**<  u wall */
    double dt = 1; /**<  just leave it a 1, simulation will crash and burn otherwise */
    double ibm_range = 2.1; /**< ibm_range */
    double lattice_length = 1.0;
    double k = 1;
    double mean_marker_distance = 0.75;
    kernelType_t kernel_in_use;
}simulation_parameters_t;

/**
 * Straight definition.
 * @note boundary type in surfaces is not yet fully supported
 */
typedef struct straight {
    // s = p + t*d
    point_t point; /**<  Origin of the straight */
    vector_t direction; /**<  Direction of the straight line */
    // validity of the straight
    double max_t = 0; /**<  Length of that line */
    boundaryType_t type = NO_BOUNDARY; /**< Boundary type */
}straight_t;

// eigen smart pointer
using link_pointer = Eigen::internal::pointer_based_stl_iterator<array_t>; /// link pointer directly to an array
using array_pointer = link_pointer; /// redefine to be more inline of what it actually does
#endif // MY_GENERAL_CODE_TYPES_H
