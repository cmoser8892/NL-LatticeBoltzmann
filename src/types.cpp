#include "types.h"

matrix_t velocity_set = {{0,1,0,-1,0 ,1,-1,-1, 1},
                         {0,0,1,0 ,-1,1, 1,-1,-1}};

matrix_t weights = {{4.0/9,1.0/9,1.0/9,1.0/9,1.0/9,1.0/36,1.0/36,1.0/36,1.0/36}};

matrix_t cardinal_directions = {{0,1,0,-1},
                                {1,0,-1,0}};

matrix_t major_directions = {{1,0,-1, 0,1,-1,-1, 1},
                             {0,1, 0,-1,1, 1,-1,-1}};