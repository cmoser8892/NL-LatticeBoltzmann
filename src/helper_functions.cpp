
#include "helper_functions.h"

bool compare_arrays(array_t a1, array_t a2) {
    if(a1.size() != a2.size()){
        return false;
    }
    else {
        for(int i = 0; i < a1.size(); ++i) {
            if(a1(i) != a2(i)) {
                return false;
            }
        }
    }
    return true;
}