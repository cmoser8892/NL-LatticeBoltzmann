/**
 * 00 is general testing ground, quick and dirty
 * when i dont want it in tests
 */


auto divide(int dividend, int divisor) {
    struct result {int quotient; int remainder;};
    return result {dividend / divisor, dividend % divisor};
}

#include <iostream>
#include "functions.h"

int main() {
    using namespace std;

    auto result = divide(14, 3);

    cout << result.quotient << ',' << result.remainder << endl;

    // or

    auto [quotient, remainder] = divide(14, 3);

    cout << quotient << ',' << remainder << endl;
}