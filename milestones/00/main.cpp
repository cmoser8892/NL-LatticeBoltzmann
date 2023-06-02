/**
 * 00 is general testing ground, quick and dirty
 * when i dont want it in tests
 */

/*
 * vallgrind call
valgrind --tool=callgrind --dump-instr=yes (p)
 */

#include <iostream>
#include "helper_functions.h"

int main() {
    solve_truncation_force_symbols(7);
}