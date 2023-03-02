#ifndef NL_LATTICEBOLTZMANN_FUNCTIONS_H
#define NL_LATTICEBOLTZMANN_FUNCTIONS_H

#include "types.h"
#include "node.h"

// functional methods
array_t equilibrium(node* node);
void collision(node* node,double relaxation);
void macro(node * node);
double bb_switch_channel(int from_channel, double uw);
int switch_link_dimensions(int link_channel);

// writeing methods
void write_rho(node* node, flowfield_t* rho);
void write_ux(node* node, flowfield_t* uy); // writes velocity based on the place in the field
void write_uy(node* node,flowfield_t* ux);
// debug method
void debug_node(node* node, bool printing);

// watchdog for rho
class rhoWatchdog {
    /// just prints out std error msgs
  private:
    flowfield_t rho;
    double sensitivity = 0.1; // displacement to the previous value
  public:
    rhoWatchdog(double s,point_t size);
    bool check(node* n,int step);
};
#endif // NL_LATTICEBOLTZMANN_FUNCTIONS_H