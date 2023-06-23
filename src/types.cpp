#include "types.h"
// common space for globals
// d2q9 velocity set
matrix_t velocity_set = {{0,1,0,-1,0 ,1,-1,-1, 1},
                         {0,0,1,0 ,-1,1, 1,-1,-1}};
// d2q9 weights
matrix_t weights = {{4.0/9,1.0/9,1.0/9,1.0/9,1.0/9,1.0/36,1.0/36,1.0/36,1.0/36}};
// cardinal directions 2d
matrix_t cardinal_directions = {{0,1,0,-1},
                                {1,0,-1,0}};