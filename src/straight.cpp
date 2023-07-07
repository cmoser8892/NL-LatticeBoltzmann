#include "straight.h"
#include "functions.h"
#include <fstream>
#include <iostream>


// private
/**
 * Calculates the mass center of all the boundary points.
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
 * Detects and moves away the mass center if it is too close to a boundary.
 */
void straightGenerator::detect_boundary_proximity_main_mass_center() {
    // checks if near a surface makes sure it is a right angle
    point_t new_mass_center = mass_center;
    for(auto surface: surfaces) {
        vector_t base = surface->point - mass_center;
        vector_t an = surface->direction;
        vector_t normal = {an.y(), - an.x()};
        double distance = base.dot(normal)/normal.norm();
        // we make use of the kernel 1 function here other kernels can be used too i guess
        double factor = kernel_1(distance,CRITICAL_RANGE);
        // we need the vector to actually move
        double mover = sqrt(base.norm()*base.norm() - distance*distance);
        point_t move_point = surface->point + an*mover;
        vector_t mover_vector = mass_center - move_point;
        // std::cout << mover_vector << std::endl;
        new_mass_center += mover_vector*factor;
        // std::cout << new_mass_center << std::endl;
    }
    mass_center = new_mass_center;
}

/**
 * Generates keys holdings for all the individual boundary-structures a boundary has.
 */
void straightGenerator::calculate_keys() {
    // fills the hash class that mirrors the boundary structures
    for(auto bs: points->boundary_structures) {
        auto pkh = new pointKeyHash;
        pkhv.push_back(pkh);
        for(auto bp : bs->boundary_points) {
            pkh->fill_key(bp->h,bp->point);
            full_pkh.fill_key(bp->h,bp->point);
        }
    }
}

/**
 * Calculates the straights between all the boundary points, only works on concave surfaces with no bumps.
 * @attention should not really matter if all the straights are saved in one place as long as different boundary structures exist
 */
void straightGenerator::calculate_all_straights() {
    // go through all the structs -> assumption closed surface
    // we just look for the nearest neighbor in the pkh in the cardinal directions
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
 * Creates all possible straights form the boundary structures.
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
        // we search in the cardinal directions for neighbors
        int found = 0;
        // we just look in the positive directions negative gets ignored :)
        // so it can happen that we do not find a partner for that point
        for(int i = 0; i < 2 ; ++i) {
            neighbor = current + point_t(cardinal_directions.col(i));
            handle_t neighbor_handle = pkh->key_translation(neighbor);
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

/**
 * Tests all the straights created in straight create and checks weather or not they are valid.
 * Combines them into longer ones if possible too.
 * @param bs of the current boundary structure
 */
void straightGenerator::straight_reduce(int bs) {
    // main function here is to reduce the number of surfaces
    // translates the temporary surface into the the used surface object
    int surface_number = 0;
    // we loop over the just created stuff
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

/**
 * We us the key hash table to find all the original boundary points to determine the length of each surface.
 * @param bs counter to the boundary structure
 */
void straightGenerator::find_surface_boundary_points(int bs) {
    // we need to check the reduced surfaces for interruptions in between
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
            if(d > 0) {
                // directional control / check
                temp /= d;
                if((temp.x() == -1) || (temp.y() == -1)) {}
                else {
                    distances.push_back(d);
                }
            };
        }
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
        current_surface->max_t = expected_distance - 1 - start_point;
        if(current_surface->max_t == 0) {
            delete current_surface;
        }
        else {
            temporary.push_back(current_surface);
        }
    }
}

/**
 * Looks for the special case of 1 long straights in the surface and resolves problems related for more details look at the note after the function header.
 * @param bs_number
 */
void straightGenerator::look_for_bumps(int bs_number) {
    /*
     * Note on what this function does:
     * - 1 long straights can be a special case with 2 corner cases associated, this function tries to tackle all of them
     * - the simplest and non corner case is two 1 opposing each other with lots of fluid between, will get ignored
     * - second we can have a various amount of 1 long straights stacked, we keep the top one and partition the bottom straight
     * - third, we get an extra 1 long straight at a corner inside the fluid, this is the special corner case that takes up most of the function to handle
     */
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
                // exclusion criteria
                if(partner == nullptr) {
                    continue;
                }
                // ray is the
                straight_t ray;
                ray.point = partner->point;
                ray.direction = partner->direction;
                double t = calculate_intersection(&ray,&surface);
                if((t >= 0) && (t <= partner->max_t)) {
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
            int total = 0;
            int expected_plus = 0;
            int expected_minus = -1;
            int position_minus_max = -1;
            int position_plus_max = -1;
            // control vector which of the values is a bad apple
            std::vector<bool> delete_true_candidates(values.size(),false);
            // loop over the found values and find the top and the bottom of the bump delete
            // prepare to delete the rest
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
                        if(position_plus_max >= 0) {
                            delete_true_candidates[position_plus_max] = true;
                        }
                        ++total;
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
                        if(position_minus_max >= 0) {
                            delete_true_candidates[position_minus_max] = true;
                        }
                        ++total;
                        position_minus_max = k;
                        --expected_minus;
                    }
                    else {
                        // set to unachievable value and safe the element
                        expected_minus = 1;
                    }
                }
            }
            // corner case we just get two participates in the surface, there is a possibility that one of them gets deleted
            bool partition_allowed = true;
            if(total == 2) {
                delete_true_candidates[0] = false;
                // we now have to determine where we are on the overall surface of the boundary
                // if we got 2 neighbors we dont have the special condition (in cardinal directions)
                // but if we got 3 we must not divide the partner on the other side
                auto pkh = pkhv[bs_number];
                // unpack
                double distance;
                double t;
                handle_t h;
                straight_t* straight;
                std::tie(distance,t,h,straight) = values[0];
                // setup and short hands
                point_t current = straight->point;
                point_t neighbor;
                handle_t neighbor_handle;
                // we search in the cardinal directions for neighbors
                int found = 0;
                // we just look in the positive directions negative gets ignored :)
                // so it can happen that we do not find a partner for that point
                for(int l = 0; l < cardinal_directions.cols() ; ++l) {
                    neighbor = current + point_t(cardinal_directions.col(l));
                    neighbor_handle = pkh->key_translation(neighbor);
                    if (neighbor_handle > 0) {
                        // self test
                        ++found;
                    }
                }
                if(found > 2) {
                    partition_allowed = false;
                }
            }
            // go over the values again and delete marked straights and part to big ones if allowed
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
                    if((straight->max_t > 1) && partition_allowed) {
                        // partition the partner
                        double lower = std::floor(t);
                        double higher = std::ceil(t);
                        // create a new surface (later part)
                        auto new_part = new straight_t;
                        new_part->point = straight->point + higher*straight->direction;
                        new_part->direction = straight->direction;
                        new_part->max_t = straight->max_t - higher;
                        // reduce the reach of the first part
                        straight->max_t = lower;
                        if(new_part->max_t > 0) {
                            temporary.push_back(new_part);
                        }
                        else {
                            delete new_part;
                        }
                        // old straight now has the length 0
                        if(straight->max_t == 0) {
                            long true_position = (long) h - 1;
                            delete temporary[true_position];
                            temporary[true_position] = nullptr;
                        }
                    }
                }
                else {
                    // delete the marked straights
                    if(delete_me) {
                        long true_position = (long) h - 1;
                        delete temporary[true_position];
                        temporary[true_position] = nullptr;
                    }
                }
            }
        }
    }
}

// public
/**
 * Constructor seed the boundary points pointer.
 * @param p
 */
straightGenerator::straightGenerator(boundaryPointConstructor *p) {
    points = p;
}

/**
 * Deletes the straights vector.
 */
straightGenerator::~straightGenerator() {
    delete_vector();
    delete_keys();
}

/**
 * Calculates the mass center and all the straights.
 */
void straightGenerator::init() {
    if(points == nullptr) {
        std::cerr << "Unintended straight generator usage" << std::endl;
        return;
    }
    // calculate main mass center and all the keys to the boundary points
    calculate_mass_center();
    // calculates hash keys for all the boundary structures
    calculate_keys();
    // see weather or not moving the mass center away from the boundaries makes sense
    // detect_boundary_proximity_main_mass_center();
    // loop over the boundary structures to create the straights
    for(int i = 0; i < points->boundary_structures.size(); ++i) {
        // we 1st create all possible straights in north and west direction
        // ( if the surface is not connected this will make it fail later on)
        straight_create(i);
        // we reduce the straights used for the surfaces to the absolute minimum
        // ( and connect in the process previously interrupted surfaces)
        straight_reduce(i);
        // we can now clear the temp creation all valids are in temp valids
        temporary_creation.clear();
        // find the boundary points of the surface again
        find_surface_boundary_points(i);
        temporary_valid.clear();
        // we look for all the small bumps created by surfaces with the length 1
        look_for_bumps(i);
        std::copy_if(temporary.begin(),
                     temporary.end(),
                     std::back_inserter(surfaces),
                     [](straight_t* s){return s != nullptr;});
        // clear temp valid too objects got added to surfaces vector
        temporary.clear();
    }
    // old legacy method
    // calculate_all_straights();
    // detect_boundary_proximity_main_mass_center();
}

/**
 * Calculates how many intersections there are between the node point and the mass center.
 * @param node_point
 * @return number of intersections
 */
int straightGenerator::calculate_intersections(const point_t node_point, point_t* individual_mc) {
    // surface based algorithm to calculate intersections
    /*
     * 3 passes have to be made to calculate to calcuate a valid intersection
     *  1 does the straight hit the surface in the area between the two points that define it
     *  2 how does the straight hit the surface (posetive or negative we only care about posetiv
     *  3 have we already hit an edgepoint
     */
    // 0 pass not a boundary point or point on the surface
    int number_of_intersections = 0;
    // determine straight to the mass center
    straight_t straight;
    straight.point = node_point; // => r
    // check if mass center and shift if yes
    if(straight.point == *individual_mc) {
        // we do a little shift out of the mass-center
        // any direction should work
        std::cerr << "Node-point is the mass-center, algorithm potentially broken" << std::endl;
        std::cerr << "Check correct node-size or amount of wet node neighbours" << std::endl;
        std::cerr << "There is just a heuristic fix to it in place" << std::endl;
        straight.point.x() += 0.1;
        straight.point.y() += 0.1;
    }
    straight.direction =  *individual_mc - straight.point;
    // setup already found
    std::vector<point_t> already_found;
    already_found.clear();
    // go over the surfaces
    for(auto surface : surfaces) {
        // in the first pass the surface is actually the ray we are using
        // we want the intersection to be between 0 and max value
        straight_t first_pass_surface;
        first_pass_surface.point = straight.point;
        first_pass_surface.direction = {straight.direction.y(), -straight.direction.x()};
        double t = calculate_intersection(surface, &first_pass_surface);
        // 1 pass
        if((t >= 0) && (t <= surface->max_t)) {
            // check if direction of the finding is positive in
            // the direction of the vector outgoing from the center
            straight_t second_pass_surface;
            second_pass_surface.point = surface->point;
            second_pass_surface.direction = {surface->direction.y(),-surface->direction.x()};
            double s = calculate_intersection(&straight,&second_pass_surface);
            // 2 pass
            if(s > 0) {
                // compare and check if we already got the same point previously
                // in case we intersect a node point directly ( no double counting )
                point_t point = straight.point + s*straight.direction;
                bool add = true;
                // 3 pass
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
            // exception the point is on the surface -> ignore/outside
            if(s == 0) {
                return 0;
            }
        }
    }
    return number_of_intersections;
}

/**
 * Supposed to be a more stable alternative for the calculate intersection function.
 * @param point
 * @return the number of intersections detected
 */
int straightGenerator::calculate_intersections_redundant(nodePoint_t *point) {
    // calculates with intersections with 3 additional mass centers in a "circle" around the original mc
    int intersection_count = 0;
    matrix_t movers = { {2,-1,-1},
                        {0, 2,-2}};
    point_t individual_mc = mass_center;
    double mover_distance = 0;
    std::vector<int> intersections;
    // calculate the intersections 3 times redundant
    for(int i = 0; i < 1; ++i) {
        individual_mc += mover_distance * (vector_t)movers.col(i);
        intersections.push_back(calculate_intersections(point->position,&individual_mc));
    }
    std::vector<int> intersections_copy = intersections;
    // we can sort and search for unique values i guess
    std::sort(intersections.begin(),intersections.end());
    auto last = std::unique(intersections.begin(),intersections.end());
    intersections.erase(last,intersections.end());
    // we check the size of the thing now
    if(intersections.size() == 1) {
        intersection_count = intersections[0];
    }
    if(intersections.size() == 2) {
        int first_number_count = std::count(intersections_copy.begin(),
                                            intersections_copy.end(),
                                            intersections[0]);
        int second_number_count = std::count(intersections_copy.begin(),
                                             intersections_copy.end(),
                                             intersections[1]);
        if(first_number_count > second_number_count) {
            intersection_count = intersections[0];
        }
        else {
            intersection_count = intersections[1];
        }
    }
    if(intersections.size() == 3) {
        // we use the original mass center here for our output
        intersection_count = calculate_intersections(point->position,&mass_center);
    }
    return intersection_count;
}

/**
 * This method is a good idea but not really implemented well.
 * @param point
 * @return
 */
bool straightGenerator::calculate_intersections_star_node_point(nodePoint_t *point) {
    // calculates intersections based on a star around the node point, no mc is needed ( we kinda cheat thou )
    // the typical roles of mc and node point are switched in this method
    // we reuse the velocity matrix for this
    bool return_value = false;
    std::vector<bool> tests;
    array_t self = point->position;
    point_t self_point = self;
    matrix_t d = {{1,0,-1,0,1,-1,-1,1},
                  {0,1,0,-1,1,1,-1,-1}};
    matrix_t e  ={{2,1,-1,-2,-2,-1,1,2},
                  {1,2,2,1,-1,-2,-2,-1}};
    matrix_t combined = d + e;
    // check in the directions of the velocity field
    for(int i = 1; i < combined.cols();++i) {
         point_t current = self + combined.col(i);
        // remember postion and mc are swapped here
        int counter = calculate_intersections(current,&self_point);
        if(counter != 0) {
            tests.push_back((counter%2) == 0);
        }
    }
    // count the numbers of true
    int passes = std::count(tests.begin(),tests.end(),true);
    if(passes > tests.size()/2) {
        return_value = true;
    }
    return return_value;
}

/**
 * Number of intersections modulo 2, if it can be divided by two and noting remains the node is inside
 * @param point
 * @return false if inside true if outside
 */
bool straightGenerator::node_inside_simple(nodePoint_t *point) {
    // even out; odd in
    // uses a surface representation to calculate weather nodes are inside or outside
    // int value = calculate_intersections(point->position, &mass_center);
    int value = calculate_intersections_redundant(point);
    return ((value%2) == 0);
}

/**
 * Node inside variant.
 * @attention Dont use, buggy not better then the current one used. Good idea though.
 * @param point
 * @return
 */
bool straightGenerator::node_inside_star(nodePoint_t *point) {
    /*
     * this method is truly crap, but the star idea is a good one
     * just not sure on the overall implementation more of a look and see feel
     */
    return calculate_intersections_star_node_point(point);
}

/**
 * Deletes the vector infos
 */
void straightGenerator::delete_vector() {
    for (auto s: surfaces) {
        delete s;
    }
}

/**
 * Delete the key classes.
 */
void straightGenerator::delete_keys() {
    for(auto k : pkhv) {
        delete k;
    }
}

/**
 * Writes out the points where the surface begins and ends xx yy.
 */
void straightGenerator::write_out_surface() {
    // function to write out the surface for boundary representation in pyplot
    std::ofstream out;
    out.open("xy_surfaces");
    if(out.is_open()) {
        // go over the surface
        for(auto s : surfaces) {
            // write out beginning and end of the surface
            point_t begin = s->point;
            point_t end = begin + s->direction*s->max_t;
            // write the info into the file
            // both x first then the ys
            out << begin.x() << "  "
                << end.x()   << "  "
                << begin.y() << "  "
                << end.y()   << "  "
                << std::endl;
        }
    }
    else {
        throw std::runtime_error("Open file has failed");
    }
}

/**
 * Adds a surface to the surface description.
 * @attention You should KNOW what you are doing here, no post checks exist for a surface
 * @param start
 * @param direction
 * @param length
 */
void straightGenerator::add_surface(straight_t s) {
    auto temp = new straight_t;
    temp->point = s.point;
    temp->direction = s.direction;
    temp->max_t = s.max_t;
    surfaces.push_back(temp);
}

/**
 * Calculates the mass_center of a given surface structure (if we add them via add_surface).
 */
void straightGenerator::surface_mass_center() {
    double point_number = 0;
    for(auto s : surfaces) {
        point_t current = s->point;
        for(int i = 0; i< s->max_t; ++i) {
            current += i*s->direction;
            mass_center += current;
            ++point_number;
        }
    }
    mass_center /= point_number;
}

/**
 * Calculates the intersection between a ray and a surface, theory is explained in the function body.
 * @param ray
 * @param surface
 * @return distance t
 */
double calculate_intersection(straight_t * ray, straight_t * surface) {
    /*
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
    point_t r = surface->point;
    vector_t n = surface->direction;
    // calculate function
    // note return nan if n orthogonal to d
    double t = ((r -o).dot(n))/(n.dot(d));
    return t;
}

/**
 * Compare function to sort tuples of values encountered in look for bumps based on distance.
 * @param a
 * @param b
 * @return
 */
bool compare_bumps_sort(const std::tuple<double,double,handle_t,straight_t*> &a,
                        const std::tuple<double,double,handle_t,straight_t*> &b) {
    // shorthand to the distance i guess
    double distance_a = std::get<0>(a);
    double distance_b = std::get<0>(b);
    // compare absolute distance values
    return std::abs(distance_a) < std::abs(distance_b);
}

/**
 * Compares the candidate with its partner returns true if the x or y value of the partner is smaller.
 * @param candidate
 * @param partner
 * @return
 */
bool straight_better_candidate_test(straight_t *candidate, straight *partner) {
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
