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
        for(int i = 0; i < bs->boundary_points.size(); ++i) {
            // search in the cardinal directions NESW
            point_t current = bs->boundary_points[array_position]->point;
            point_t candidate;
            handle_t candidate_handle = 0;
            for(int i = 0; i < 4; ++i) {
                candidate = current + point_t(cardinal_directions.col(i));
                // std::cout << candidate << std::endl;
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
            // rotate the vector by 90 degrees forward (doesnt really matter which direction)
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
    calculate_all_straights();
}

void straightGenerator::init_test() {
    // just for testing out the new init structure
    // will get added in the boundary
    // 1st step generate all connections in immediate neighborhood of a point

}

/**
 * @fn int straightGenerator::calculate_intersections(nodePoint_t* node_point)
 * @brief calculates how many intersections there are between the node point and the mass center
 * @param node_point
 * @return number of intersections
 */
int straightGenerator::calculate_intersections(nodePoint_t* node_point) {
    /// surface based algorithm to calcualte intersections
    /**
     * 3 passes have to be made to calculate to calcuate a valid intersection
     *  1 does the straight hit the surface in the area between the two points that define it
     *  2 how does the straight hit the surface (posetive or negative we only care about posetiv
     *  3 have we already hit an edgepoint
     */
    int number_of_intersections = 0;
    // check if actually the boundary point
    for(auto surf : surfaces) {
        // check if we are a surface point described (aka a boundary point and do a hard break
        if(surf->point == point_t(node_point->position)) {
            return 0;
        }
    }
    // corner case mass center lays on point
    point_t p = node_point->position;
    // determine straight to the mass center
    straight_t straight;
    straight.point = node_point->position; // => r
    // check if mass center
    if(straight.point == mass_center) {
        // we do a little shift out of the mass-center
        // any direction should work
        std::cerr << "Node-point is the mass-center, algorithm potentially broken" << std::endl;
        std::cerr << "Check correct node-size or amount of wet node neighbours" << std::endl;
        straight.point.x() += 0.1;
        straight.point.y() += 0.1;
    }
    straight.direction =  mass_center - straight.point;
    vector_t normal = {straight.direction.y(), -straight.direction.x()}; // => n
    // go through the surface and take a look
    std::vector<point_t> already_found;
    already_found.clear();
    // actual test
    for(auto surf : surfaces) {
        // t = ((r - o)·n)/(n·d)
        // surf->point => o
        // surf->direction => d
        // std::cout << straight.point << std::endl;
        double t = ((straight.point - surf->point).dot(normal))/
                   (normal.dot(surf->direction));
        if((t >= 0.0) && (t <= 1.0)) {
            // check if direction of the finding is posetiv in the direction of the vector
            vector_t surface_normal = {surf->direction.y(),-surf->direction.x()};
            double s = ((surf->point-straight.point).dot(surface_normal))/(surface_normal.dot(straight.direction));
            if(s >= 0) {
                point_t point = straight.point + s*straight.direction;
                bool add = true;
                // might be a bad idea no just check s and not o +sd
                // but it should be alright
                for (auto ps : already_found) {
                    if(ps == point) {
                        add = false;
                    }
                }
                already_found.push_back(point);
                if(add) {
                    number_of_intersections++;
                }
            }
        }
    }
    if(0) {
        std::cout << "Result" << std::endl;
        std::cout << node_point->position.x() << " " << node_point->position.y() << std::endl;
        std::cout << number_of_intersections << std::endl;
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

