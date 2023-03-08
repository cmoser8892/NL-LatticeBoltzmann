#include <iostream>
#include "simulation.h"

/**
*
 * todos write more docu, do this lazylly at home thou
 *
 * Todos/Findings:
    looked through the big callgrids, its quite insane how much collision takes
    calculate seems superior than the lookup thou, prob comes down to cpu stalls, should also do a cache hit analysis prob to shed some light on that
    breakdown of rho seems to be systematic and not related to my code could reproduce it also in the older python code, so i will just accept it, my inti thinks its something with momentums in the moving bb
*/

int main(int argc, char *argv[]) {
    std::cout << "hi" << std::endl;
}
