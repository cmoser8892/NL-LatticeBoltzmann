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
    // todo fix errors related to mass center placement in the check inside algorithm :)
    also todo  make 3 mcs for a more stable algorithm
    how do i test sth to be convex?!
    write one more test related to the construction oposing bumps should not be a prob thou
    mass_center /=3;
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
        // we just look in the positive directions negative gets ignored :)
        // so it can happen that we do not find a partner for that point
        for(int i = 0; i < 2 ; ++i) {
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
    }
}

void straightGenerator::straight_self_test(int bs) {
    // translates the temporary surface into the the used surface object
    int surface_number = 0;
    for(int i = 0; i<temporary_creation.size(); ++i) {
        auto candidate = temporary_creation[i];
        if(candidate == nullptr) {
            continue;
        }
        bool better_candidate = false;
        // shorthands
        point_t current = candidate->point;
        vector_t direction = candidate->direction;
        // set up a straight line through the middle of the candidate normal to the direction
        straight_t ray;
        ray.point = current + 0.5 * direction;
        ray.direction = {direction.y(), -direction.x()};
        // we do an intersection test with all the other surfaces in the temporary object
        for(int j = 0; j < temporary_creation.size();++j) {
            auto partner = temporary_creation[j];
            // ignore self
            if(partner == candidate) {
                continue;
            }
            // ignore deleted
            if(partner == nullptr) {
                continue;
            }
            // set up the surface
            straight_t surface;
            surface.point = partner->point;
            surface.direction = {partner->direction.y(), -partner->direction.x()};
            // do an intersection test
            // to be more sure of the result we round
            double t = std::round(calculate_intersection(&ray,&surface));
            if(t == 0) {
                // found an equal surface delete that one
                // if a point is further starts before the positive direction point elect that one and
                // stop this loop to delete the candidate
                if(!straight_better_candidate_test(candidate,partner)) {
                    delete partner;
                    temporary_creation[j] = nullptr;
                }
                else {
                    better_candidate = true;
                    break;
                }
            }
        }
        // add to valid if nothing better is found
        if(!better_candidate) {
            temporary_valid.push_back(candidate);
        }
        else {
            delete candidate;
            temporary_creation[i] = nullptr;
        }
    }
}

bool straightGenerator::straight_closer_test(int bs ,straight_t* self, straight_t* partner) {
    bool return_value = false;
    bool distance_test = false;
    // make sure we are in the boarders of the other guy
    // the ray here is given by the parter
    straight_t surface_self;
    surface_self.point = self->point + 0.5*self->direction;
    surface_self.direction = {self->direction.y(), -self->direction.x()};
    double t = calculate_intersection(partner, &surface_self);
    if((t >= 0.0) && (t <= 1.0)) {
        distance_test = true;
    }
    // do the distance test
    if(distance_test) {
        // setup for who is more distant to the mc
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
    }
    return return_value;
}

bool straightGenerator::straight_better_candidate_test(straight_t *candidate, straight *partner) {
    bool return_value = false;
    // if either x or y value of the origin is lower prob better
    if(partner->point.x() < candidate->point.x() ) {
        return_value = true;
    }
    if(partner->point.y() < candidate->point.y()) {
        return_value = true;
    }
    return return_value;
}

void straightGenerator::find_surface_boundary_points(int bs) {
    auto bps = points->boundary_structures[bs]->boundary_points;
    for (auto s : temporary_valid) {
        std::vector<point_t*> on_surface_points;
        std::vector<double> distances;
        vector_t direction = {s->direction.y(),-s->direction.x()};
        for( auto bp : bps) {
            // construct
            vector_t check = bp->point - s->point;
            if(std::round(direction.dot(check)) == 0) {
                // orthogonal so on the surface
                on_surface_points.push_back(&bp->point);
            }
        }
        // we found all the boundary points on the structure
        // now we determine the t values
        // first we determine all the distances form the origin of the surface
        for(auto surface_points : on_surface_points) {
            vector_t temp = (*surface_points - s->point);
            double d = temp.norm();
            // std::cout << d << std::endl;
            if(d > 0) {
                // directional control / check
                temp /= d;
                if((temp.x() == -1) || (temp.y() == -1)) {
                    /// todo negate this expression to be a oneliner
                }
                else {
                    distances.push_back(d);
                }
            };
        }
        // std::cout << std::endl;
        // then we make sure to have them sorted
        std::sort(distances.begin(),distances.end());
        // now we use that to determine t values and check weather or not we have to partition the thing
        // control variables
        double expected_distance = 1;
        double start_point = 0;
        int iteration_counter = 0;
        auto current_surface = s;
        point_t current_surface_start_point = s->point;
        int i = 0;
        // do while loop for easier control
        while(i < distances.size()){
            // setup d
            auto d = distances[i];
            // ignore negative distances
            if(d > 0) {
                //
                if(d != expected_distance) {
                    // first entry is more away than expected
                    if(iteration_counter == 0) {
                        // error ?!
                    }
                    // set t values
                    current_surface->min_t = 0;
                    current_surface->max_t = expected_distance - 1 - start_point;
                    // push into the temporary field
                    temporary.push_back(current_surface);
                    // create and setup a new surface
                    auto temp = new straight_t;
                    temp->point = current_surface_start_point + d*current_surface->direction;
                    temp->direction = current_surface->direction;
                    current_surface = temp;
                    // reset control variables
                    expected_distance = d + 1;
                    iteration_counter = 0;
                    start_point = d;
                }
                else {
                    ++expected_distance;
                    ++iteration_counter;
                }
            }
            // iterator up
            ++i;
        }
        // set values
        current_surface->min_t = 0;
        current_surface->max_t = expected_distance - 1 - start_point;
        if((current_surface->min_t == 0) && (current_surface->max_t == 0)) {
            delete current_surface;
        }
        else {
            temporary.push_back(current_surface);
        }
    }
}

void straightGenerator::straight_set_t_values(int bs_number) {
    // set up the correct structures
    auto pkh = pkhv[bs_number];
    // we go through temp_valid and look in minus and plus directions to set the t values
    for(auto s : temporary_valid) {
        // set min_t
        s->min_t = go_through_vector(bs_number,s,-1);
        s->max_t = go_through_vector(bs_number,s,1);
        // if the surface is good one should not be 0
        if((s->min_t == 0) && (s->max_t == 0)) {
            throw std::runtime_error("bad surface");
        }
        surfaces.push_back(s);
    }
}

double straightGenerator::go_through_vector(int bs_number, straight_t *self, int plus_minus) {
    double return_value = 0;
    // set up the correct structures
    auto pkh = pkhv[bs_number];
    long max_iterations = points->total_boundary_nodes();
    // shorthands
    point_t origin = self->point;
    point_t current = origin;
    point_t direction = self->direction;
    int previous = 0;
    for(int i = 0; i < max_iterations; ++i) {
        current = origin + direction*plus_minus*i;
        if(pkh->key_translation(current)) {
            previous = i;
        }
        else {
            return_value = plus_minus*previous;
            break;
        }
    }
    return return_value;
}

void straightGenerator::look_for_bumps(int bs) {
    // test out straights with length 1
    for(int i = 0; i < temporary.size();++i) {
        // anti nullptr work
        auto self = temporary[i];
        if(self == nullptr) {
            continue;
        }
        // only do something if we are the length 1
        if(self->max_t == 1) {
            // we do an intersection test with all the others
            // if we intersect between 0 and max we have to decide weather this one is necessary
            // we also have to check the distance between them
            // all the distances get recorded too.
            // construct the surface, just move the origin to the middle
            straight_t surface;
            surface.point = self->point + 0.5 * self->direction;
            surface.direction = self->direction;
            // setup vectors
            std::vector<std::tuple<double,double,handle_t, straight_t*>> values;
            // search for potential parters
            handle_t handle = 0;
            for(auto partner : temporary) {
                ++handle;
                // exclusion craterias
                if(partner == nullptr) {
                    continue;
                }
                // ray is the
                straight_t ray;
                ray.point = partner->point;
                ray.direction = partner->direction;
                double t = calculate_intersection(&ray,&surface);
                if((t >= partner->min_t) && (t <= partner->max_t)) {
                    // why does eigen not support 2d cross products?!
                    // set up everything to work correctlly
                    straight_t partner_surface;
                    partner_surface.point = partner->point;
                    partner_surface.direction = {partner->direction.y(), - partner->direction.x()};
                    straight_t self_ray;
                    self_ray.point = self->point;
                    self_ray.direction = {self->direction.y(), -self->direction.x()};
                    double s = calculate_intersection(&self_ray,&partner_surface);
                    // add to lists
                    values.push_back(std::make_tuple(s,t,handle,partner));
                }
            }
            // if values has just the size 2 we can skip the next steps
            if(values.size() == 2) {
                continue;
            }
            // sort via the abs values of the distance
            std::sort(values.begin(),values.end(), compare_bumps_sort);
            // go through list of partners
            int expected_plus = 0;
            int expected_minus = -1;
            int position_minus_max = -1;
            int position_plus_max = -1;
            // control vector which of the values is a bad apple
            std::vector<bool> delete_true_candidates(values.size(),false);
            if(values.size() > 2) {
                std::cout << "here" << std::endl;
            }
            // loop over the found values
            for(int k = 0; k < values.size();++k) {
                auto candy = values[k];
                // unpack
                double distance;
                double t;
                handle_t h;
                straight_t* straight;
                std::tie(distance,t,h,straight) = candy;
                // validity test
                if(distance >= 0) {
                    // positive
                    if (distance == expected_plus) {
                        // enable delete for the previous candidate
                        if(position_plus_max >= 0)
                            delete_true_candidates[position_plus_max] = true;

                        position_plus_max = k;
                        ++expected_plus;
                    }
                    else {
                        // set to unclear values and save the element
                        expected_plus = -1;
                    }
                }
                else {
                    // negative
                    if(distance == expected_minus) {
                        // enable delete for the previous candidate
                        if(position_minus_max >= 0)
                            delete_true_candidates[position_minus_max] = true;
                        position_minus_max = k;
                        --expected_minus;
                    }
                    else {
                        // set to unachievable value and safe the element
                        expected_minus = 1;
                    }
                }
            }
            for(int k = 0; k < values.size(); ++k) {
                auto candy = values[k];
                auto delete_me = delete_true_candidates[k];
                // unpack
                double distance;
                double t;
                handle_t h;
                straight_t* straight;
                std::tie(distance,t,h,straight) = candy;
                if((position_minus_max == k) || (position_plus_max == k)) {
                    if(straight->max_t > 1) {
                        // partition the partner
                        double lower = std::floor(t);
                        double higher = std::ceil(t);
                        // create a new surface (later part)
                        auto new_part = new straight_t;
                        new_part->point = straight->point + higher*straight->direction;
                        new_part->direction = straight->direction;
                        new_part->min_t =  0;
                        new_part->max_t = straight->max_t - higher;
                        // reduce the reach of the first part
                        if(new_part->max_t > 0) {
                            straight->max_t = lower;
                            temporary.push_back(new_part);
                        }
                        else {
                            delete new_part;
                        }
                    }
                }
                else {
                    // todo not that simple there are cases where the opposite side
                    // is just 1 long!!!
                    if(delete_me) {
                        long true_position = (long) h - 1;
                        delete temporary[true_position];
                        temporary[true_position] = nullptr;
                    }
                }
            }
            /*
            if(0) {
                // completely independent of expected distance only necessary
                // to be a valid one
                if((straight->max_t == 1) && (valid)) {
                    // make sure who is more close to mc?!
                    delete temporary[i];
                    temporary[i] = nullptr;
                }
                else if((straight->max_t > 1) && (valid)) {
                    // partition the partner
                    double lower = std::floor(t);
                    double higher = std::ceil(t);
                    // create a new surface (later part)
                    auto new_part = new straight_t;
                    new_part->point = straight->point + higher*straight->direction;
                    new_part->direction = straight->direction;
                    new_part->min_t =  0;
                    new_part->max_t = straight->max_t - higher;
                    // reduce the reach of the first part
                    straight->max_t = lower;
                    temporary.push_back(new_part);
                }
            }
             */
        }
    }
}

void straightGenerator::straight_test_creation(int bs) {
    for(auto s : temporary) {
        if(s != nullptr) {
            surfaces.push_back(s);
        }
    }
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
    delete_keys();
}

/**
 * @fn void straightGenerator::init()
 * @brief calculates the mass center and all the straights
 */
void straightGenerator::init() {
    calculate_mass_center();
    calculate_keys();
    init_test();
    // calculate_all_straights();
}

void straightGenerator::init_test() {
    // just for testing out the new init structure
    // will get added in the boundary
    // 1st step generate all connections in immediate neighborhood of a point
    for(int i = 0; i < points->boundary_structures.size(); ++i) {
        straight_create(i);
        straight_self_test(i);
        // we can now clear the temp creation all valids are in temp valids
        temporary_creation.clear();
        find_surface_boundary_points(i);
        temporary_valid.clear();
        // todo look for bumps is a bad function breaks more than just one test
        look_for_bumps(i);
        straight_test_creation(i);
        // clear temp valid too objects got added to surfaces vector
        temporary.clear();
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
    for(auto bs: points->boundary_structures) {
        for(auto bp : bs->boundary_points) {
            // check if we are a surface point described (aka a boundary point and do a hard break
            if(bp->point == point_t(node_point->position)) {
                return 0;
            }
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
        if((t >= surface->min_t) && (t <= surface->max_t)) {
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

void straightGenerator::delete_keys() {
    for(auto k : pkhv) {
        delete k;
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
     *   where o and d are the direction of a ray
     *   and r and n describe a surface point and a normal to describe a surface
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

bool compare_bumps_sort(const std::tuple<double,double,handle_t,straight_t*> &a,
                        const std::tuple<double,double,handle_t,straight_t*> &b) {
    // shorthand to the distance i guess
    double distance_a = std::get<0>(a);
    double distance_b = std::get<0>(b);
    // compare absolute distance values
    return std::abs(distance_a) < std::abs(distance_b);
}