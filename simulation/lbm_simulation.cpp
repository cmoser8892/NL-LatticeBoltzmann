# include "lbm_simulation.h"

std::tuple<double, double, double> ibmSimulation::calcualte_macro(array_t *a, array_t *previous) {
    return std::tuple<double, double, double>();
}
void ibmSimulation::streaming(array_t *a, std::vector<link_pointer> *list) {
}
void ibmSimulation::collision(array_t *a, double rho) {
}
void ibmSimulation::forcing_term() {
}
void ibmSimulation::lbm_interpolate_forward_compute() {
}
ibmSimulation::ibmSimulation(nodeGenerator *g, goaForce *f, markerIBM *m, vector_t size) {
}
ibmSimulation::~ibmSimulation() {
}
void ibmSimulation::set_simulation_parameters(simulation_parameters_t t) {
}
void ibmSimulation::init() {
}
void ibmSimulation::run(int current_step) {
}
void ibmSimulation::get_data(bool write_to_file) {
}
void ibmSimulation::delete_nodes() {
}
