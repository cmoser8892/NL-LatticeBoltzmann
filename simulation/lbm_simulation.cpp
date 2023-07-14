# include "lbm_simulation.h"

std::tuple<double, double, double> ibmSimulation::calculate_macro(array_t *a, array_t *previous) {
    return std::tuple<double, double, double>();
}
void ibmSimulation::streaming(array_t *a, std::vector<link_pointer> *list) {
    // just the sim
    for(int i = 1; i < CHANNELS; ++i) {
        // pointer magic
        auto origin = a->begin() + offset_node + i;
        (list->operator[](i-1) + offset_sim).operator*() = origin.operator*();
    }
}
void ibmSimulation::collision(array_t *a, double rho, double ux, double uy) {
    // undrosed collision term
    int o = offset_node;
    double relaxation = 1/parameters.relaxation;
    auto p = a->begin() + o;
    (p + 0).operator*() -= relaxation * ((p + 0).operator*() - weights.col(0).x()*rho*(1- 1.5*(ux*ux +uy*uy)));
    (p + 1).operator*() -= relaxation * ((p + 1).operator*() - weights.col(1).x()*rho*(1+ 3*ux+ 4.5*ux*ux- 1.5*(ux*ux +uy*uy)));
    (p + 2).operator*() -= relaxation * ((p + 2).operator*() - weights.col(2).x()*rho*(1+ 3*uy+ 4.5*uy*uy- 1.5*(ux*ux +uy*uy)));
    (p + 3).operator*() -= relaxation * ((p + 3).operator*() - weights.col(3).x()*rho*(1- 3*ux+ 4.5*ux*ux- 1.5*(ux*ux +uy*uy)));
    (p + 4).operator*() -= relaxation * ((p + 4).operator*() - weights.col(4).x()*rho*(1- 3*uy+ 4.5*uy*uy- 1.5*(ux*ux +uy*uy)));
    (p + 5).operator*() -= relaxation * ((p + 5).operator*() - weights.col(5).x()*rho*(1+ 3*ux+ 3*uy+ 9*ux*uy+ 3*(ux*ux +uy*uy)));
    (p + 6).operator*() -= relaxation * ((p + 6).operator*() - weights.col(6).x()*rho*(1- 3*ux+ 3*uy- 9*ux*uy+ 3*(ux*ux +uy*uy)));
    (p + 7).operator*() -= relaxation * ((p + 7).operator*() - weights.col(7).x()*rho*(1- 3*ux- 3*uy+ 9*ux*uy+ 3*(ux*ux +uy*uy)));
    (p + 8).operator*() -= relaxation * ((p + 8).operator*() - weights.col(8).x()*rho*(1+ 3*ux- 3*uy- 9*ux*uy+ 3*(ux*ux +uy*uy)));
}
void ibmSimulation::forcing_term() {
}
void ibmSimulation::lbm_interpolate_forward_compute() {
}
ibmSimulation::ibmSimulation(nodeGenerator *g, goaForce *f, markerIBM *m, vector_t size) {
}
ibmSimulation::~ibmSimulation() {
    delete_containers();
}
void ibmSimulation::set_simulation_parameters(simulation_parameters_t t) {
}
void ibmSimulation::init() {
}
void ibmSimulation::run(int current_step) {
}
void ibmSimulation::get_data(bool write_to_file) {
}
void ibmSimulation::delete_containers() {
    for(auto n : nodes) {
        delete n;
    }
    for(auto m : markers) {
        delete m;
    }
    nodes.clear();
    markers.clear();
}
