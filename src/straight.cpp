#include "straight.h"
#include <iostream>

/// private
/**
 * @fn void straightGenerator::calculate_mass_center()
 * @brief calculates the mass center of all the boundary points
 */
void straightGenerator::calculate_mass_center() {
    // add all the points up, divide by the number of points
    // as an alternative the individual mcs could be added up
    mass_center.setZero();
    for(auto bs : points->boundary_structures) {
        // set up the individual mass center
        point_t imc;
        imc.setZero();
        // add up
        for(auto b : bs->boundary_points) {
            imc += b->point;
            mass_center += b->point;
        }
        // divide by the total size
        imc /= double(bs->boundary_points.size());
        // save the individual mc of the boundary_structures
        individual_mass_centers.push_back(imc);
    }
    // div by the total size
    mass_center /= double(points->total_boundary_nodes());
}

/**
 * @fn void straightGenerator::calculate_keys()
 * @brief generates keys holdings for all the individual boundary-structures a boundary has
 */
void straightGenerator::calculate_keys() {
    for(auto bs: points->boundary_structures) {
        auto pkh = new pointKeyHash;
        pkhv.push_back(pkh);
        for(auto bp : bs->boundary_points) {
            pkh->fill_key(bp->h,bp->point);
        }
    }
}

/**
 * @fn void straightGenerator::calculate_all_straights()
 * @brief calculates the straights between all the boundary points
 * @attention should not really matter if all the straights are saved in one place as long as different boundary structures exist
 */
void straightGenerator::calculate_all_straights() {
    // go through all the structs -> assumption closed surface
    for(int i = 0; i< points->boundary_structures.size(); ++i) {
        auto bs = points->boundary_structures[i];
        auto pkh = pkhv[i];
        std::vector<handle_t> blacklist;
        // first point
        handle_t array_position = 0;
        for(int j = 0; j < bs->boundary_points.size(); ++j) {
            // search in the cardinal directions NESW
            point_t current = bs->boundary_points[array_position]->point;
            point_t candidate;
            handle_t candidate_handle = 0;
            for(int k = 0; k < cardinal_directions.cols()  ; ++k) {
                candidate = current + point_t(cardinal_directions.col(i));
                // linear search over all keys lol
                candidate_handle = pkh->key_translation(candidate);
                // if not 0 this is valid handle
                if(candidate_handle > 0) {
                    // search form the back (more likely the last two elements)
                    bool already_visited = false;
                    auto iter = blacklist.rbegin();
                    while(iter != blacklist.rend()) {
                        if(iter.operator*() == candidate_handle) {
                            // bad egg found
                            already_visited = true;
                            break;
                        }
                        iter++;
                    }
                    // break general search in cardinal directions viable candidate found
                    if(!already_visited) {
                        break;
                    }
                }
            }
            if(candidate_handle == 0) {
                throw std::runtime_error("construction of the containing surface failed!");
            }
            // add the surface
            auto s  = new straight_t;
            s->point = current;
            // calculate the vector
            vector_t next = candidate - current;
            s->direction = next;
            // put the point in the middle
            surfaces.push_back(s);
            // add to blacklist
            blacklist.push_back(candidate_handle);
            // check for overflow!
            array_position = candidate_handle-1;
        }
    }
}

/**
 * @fn void straightGenerator::straight_create(int bs_number)
 * @brief creates all possible
 * @param bs_number (which boundary structure based on index)
 */
void straightGenerator::straight_create(int bs_number) {
    // creates a temporary object based for the boundary_structre
    auto bs = points->boundary_structures[bs_number];
    auto pkh = pkhv[bs_number];
    for(auto bp : bs->boundary_points) {
        // setup and short hands
        point_t current = bp->point;
        point_t neighbor;
        handle_t neighbor_handle;
        // we search in the cardinal directions for neighbors
        int found = 0;
        for(int i = 0; i < cardinal_directions.cols()  ; ++i) {
            neighbor = current + point_t(cardinal_directions.col(i));
            neighbor_handle = pkh->key_translation(neighbor);
            if(neighbor_handle > 0) {
                // we add it to the temporary object used in self test
                auto s  = new straight_t;
                s->point = current;
                s->direction = neighbor - current;
                temporary_creation.push_back(s);
                // self test
                ++found;
            }
        }
        if(found == 0) {
            throw std::runtime_error("construction of the containing surface failed!");
        }
    }
}

void straightGenerator::straight_self_test(int bs) {
    // translates the temporary surface into the the used surface object
    for(auto candidate :temporary_creation) {
        // control variable
        bool add_me = true;
        // shorthands
        point_t current = candidate->point;
        vector_t direction = candidate->direction;
        // set up a straight line through the middle of the candidate normal to the direction
        straight_t ray;
        ray.point = current + 0.5 * direction;
        ray.direction = {direction.y(), -direction.x()};
        // we do an intersection test with all the other surfaces in the temporary object
        for(auto partner : temporary_creation) {
            // ignore self
            if(partner == candidate) {
                continue;
            }
            // set up the surface
            straight_t surface;
            surface.point = partner->point;
            surface.direction = {partner->point.y(), -partner->point.x()};
            // do an intersection test
            // to be more sure of the result we round
            double t = std::round(calculate_intersection(&ray,&surface));
            if(t == 0) {
                // found an equal surface dont add
                add_me = false;
                break;
            }
            if((t == 1) ||(t == -1)) {
                // found an adjacent surface
                // we now test who is closer to the mass center and set add_me according to that
                if(straight_closer_test(bs,candidate,partner)) {
                    add_me = false;
                    break;
                }
            }
        }
        if(add_me) {
            // add to the general surface
            surfaces.push_back(candidate);
        }
        else {
            // delete not necessary
            delete candidate;
        }
    }
}

bool straightGenerator::straight_closer_test(int bs ,straight_t* self, straight_t* partner) {
    bool return_value = false;
    // setup
    point_t current_mc = individual_mass_centers[bs];
    vector_t mc_self_midpoint = (self->point + 0.5 * self->direction) - current_mc;
    vector_t mc_partner_midpoint = (partner->point + 0.5 * partner->direction) - current_mc;
    // which one is longer ?!
    if(mc_partner_midpoint.norm() > mc_self_midpoint.norm()) {
        return_value = true;
    }
    else if(mc_self_midpoint.norm() == mc_partner_midpoint.norm()) {
        throw std::runtime_error("Algorithm can not decide which part of the surface to discard");
    }
    return return_value;
}

/// public
/**
 * @fn straightGenerator::straightGenerator(boundaryPointConstructor *p)
 * @brief constructor seed the boundary points pointer
 * @param p
 */
straightGenerator::straightGenerator(boundaryPointConstructor *p) {
    points = p;
}

/**
 * @fn straightGenerator::~straightGenerator()
 * @brief deletes the straights vector
 */
straightGenerator::~straightGenerator() {
    delete_vector();
}

/**
 * @fn void straightGenerator::init()
 * @brief calculates the mass center and all the straights
 */
void straightGenerator::init() {
    calculate_mass_center();
    calculate_keys();
    init_test();
    calculate_all_straights();
}

void straightGenerator::init_test() {
    // just for testing out the new init structure
    // will get added in the boundary
    // 1st step generate all connections in immediate neighborhood of a point
    for(int i = 0; i < points->boundary_structures.size(); ++i) {
        straight_create(i);
        straight_self_test(i);
    }
}

/**
 * @fn int straightGenerator::calculate_intersections(nodePoint_t* node_point)
 * @brief calculates how many intersections there are between the node point and the mass center
 * @param node_point
 * @return number of intersections
 */
int straightGenerator::calculate_intersections(nodePoint_t* node_point) {
    /// todo why 3 passes should not 2 be enough
    /// surface based algorithm to calculate intersections
    /**
     * 3 passes have to be made to calculate to calcuate a valid intersection
     *  0 not a boundary point used in construction
     *  1 does the straight hit the surface in the area between the two points that define it
     *  2 how does the straight hit the surface (posetive or negative we only care about posetiv
     *  3 have we already hit an edgepoint
     */
    /// 0 pass not a boundary point or point on the surface
    int number_of_intersections = 0;
    // check if actually the boundary point, boundary points are excluded in the first pass
    for(auto surf : surfaces) {
        // check if we are a surface point described (aka a boundary point and do a hard break
        if(surf->point == point_t(node_point->position)) {
            return 0;
        }
    }
    // determine straight to the mass center
    straight_t straight;
    straight.point = node_point->position; // => r
    // check if mass center and shift if yes
    if(straight.point == mass_center) {
        // we do a little shift out of the mass-center
        // any direction should work
        std::cerr << "Node-point is the mass-center, algorithm potentially broken" << std::endl;
        std::cerr << "Check correct node-size or amount of wet node neighbours" << std::endl;
        straight.point.x() += 0.1;
        straight.point.y() += 0.1;
    }
    straight.direction =  mass_center - straight.point;
    // setup already found
    std::vector<point_t> already_found;
    already_found.clear();
    // go over the surfaces
    for(auto surface : surfaces) {
        // in the first pass the surface is actually the ray we are using
        // we want the intersection to be between 0 and 1
        straight_t first_pass_surface;
        first_pass_surface.point = straight.point;
        first_pass_surface.direction = {straight.direction.y(), -straight.direction.x()};
        double t = calculate_intersection(surface, &first_pass_surface);
        /// 1 pass
        if((t >= 0.0) && (t <= 1.0)) {
            // check if direction of the finding is positive in
            // the direction of the vector outgoing from the center
            straight_t second_pass_surface;
            second_pass_surface.point = surface->point;
            second_pass_surface.direction = {surface->direction.y(),-surface->direction.x()};
            double s = calculate_intersection(&straight,&second_pass_surface);
            /// 2 pass
            if(s >= 0) {
                // compare and check if we already got the same point previously
                // in case we intersect a node point directly ( no double counting )
                point_t point = straight.point + s*straight.direction;
                bool add = true;
                /// 3 pass
                for (auto&  ps : already_found) {
                    if(ps == point) {
                        add = false;
                    }
                }
                // push the point into the already found bin
                already_found.push_back(point);
                // increase the numer of intersections
                if(add) {
                    number_of_intersections++;
                }
            }
        }
    }
    return number_of_intersections;
}

/**
 * @fn bool straightGenerator::node_inside(nodePoint_t *point)
 * @brief number of intersections modulo 2, if it can be divided by two and noting remains the node is inside
 * @param point
 * @return false if inside true if outside
 */
bool straightGenerator::node_inside(nodePoint_t *point) {
    // even out; odd in
    /// uses a surface representation to calculate weather nodes are inside or outside
    int value = calculate_intersections(point);
    return ((value%2) == 0);
}

/**
 * @fn void straightGenerator::delete_vector()
 * @brief deletes the vector infos
 */
void straightGenerator::delete_vector() {
    for (auto s: surfaces) {
        delete s;
    }
}

/**
 * @fn double calculate_intersection(straight_t * ray, straight_t * surface)
 * @brief calculates the intersection between a ray and a surface, theory is explained in the function body
 * @param ray
 * @param surface
 * @return distance t
 */
double calculate_intersection(straight_t * ray, straight_t * surface) {
    /**
     * Theory:
     * - An intersection occurs if the equation f(x,y) = f(r(t)) = f(o + t*d) = 0 is satisfied
     * - All points p = (x,y) on a plane
     *   with surface normal n and offset r satisfy the equation n dot (p - r) = 0
     * - Therefore the intersection with the ray can be computed based on t:
     *   n dot (o + t*d - r) = 0 --> t = ((r - o) dot n)/(n dot d)
     * - we use t to decide on the type of intersection and weather or not it is valid in
     *   the specific use case
     */
    // shorthands
    point_t o = ray->point;
    vector_t d = ray->direction;
    // it is assumed that the surface normal is given
    // todo save the surface normals instead!
    point_t r = surface->point;
    vector_t n = surface->direction;
    // calculate function
    /// note return nan if n orthogonal to d
    double t = ((r -o).dot(n))/(n.dot(d));
    return t;
}