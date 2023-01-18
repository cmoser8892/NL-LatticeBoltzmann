#include <iostream>
#include "types.h"
using namespace std;


int main(){
  cout << "hi" << endl;
  point_t o = {3,1};
  vector_t d = {1,0};
  point_t r = {5,-1};
  vector_t n = {-1,-1};
  double t = ((r - o).dot(n)) / (d.dot(n));
  cout << t << endl;
  cout << o + d*t << endl;
}
