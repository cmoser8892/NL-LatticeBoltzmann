#include "nodeGenerator.h"
#include "helper_functions.h"
#include "simulation.h"
#include "neighborhood.h"
#include <iostream>
#include <fstream>
#include <filesystem>

/**
 * @fn nodeGenerator::nodeGenerator(boundaryPointConstructor *p)
 * @brief constructor, sets points
 * @param p
 */
nodeGenerator::nodeGenerator(boundaryPointConstructor *p) {
    points = p;
}

/**
 * @fn nodeGenerator::~nodeGenerator()
 * @brief deletes all the info about the nodes
 */
nodeGenerator::~nodeGenerator() {
    delete_node_infos();
}
/**
 * @fn void nodeGenerator::linear_generation()
 * @brief generates nodes until it hits a boundary useful to make simple 1d tests
 */
void nodeGenerator::linear_generation() {
    handle_t handle_counter = 1;
    // go throu the boundary points starting at a b point and go throu while still discovering new ones
    for(auto p : points->boundary_points) {
        point_t current = p->point + discovery_vector;
        // until a boundary points is hit create more points
        // also check weather or not inside boundaries
        while(!check_other_boundary_hit(p,current)) {
            auto n = new nodePoint_t;
            n->handle = handle_counter;
            n->position = current;
            n->type = WET;
            n->boundary = NO_BOUNDARY;
            // dont forget to increase the handle counter each time
            handle_counter++;
            node_infos.push_back(n);
            current += discovery_vector;
        }
    }
    // finally add boundary nodes
    add_boundary_nodes(&handle_counter);
}
/**
 * @fn bool nodeGenerator::check_other_boundary_hit(boundaryPoint_t* p,point_t &check_point)
 * @brief checks if the next point would be a boundary or outside
 * @param p
 * @param check_point
 * @return true or false based if boundary/not also checks if still inside bounds
 */
bool nodeGenerator::check_other_boundary_hit(boundaryPoint_t* p,point_t &check_point) {
    // returns true if hit or outside false if not
    bool return_value = false;
    point_t check = check_point.base();
    // right now just does a linear search to check if point hit or not
    for(auto point : points->boundary_points) {
        if(compare_two_points(&check,&point->point)) {
            return_value = true;
            break;
        }
    }
    // then check if still inside sim space
    point_t limit_lower = {0,0};
    point_t limit_upper = points->limits;
    if(!check_inside_limits_upper_lower(&check,&limit_lower,&limit_upper)) {
        return_value = true;
    }
    return return_value;
}


/**
 * @fn void nodeGenerator::determine_neighbors()
 * @brief hull for the neighbourhood class with its determine_neighbors call
 */
void nodeGenerator::determine_neighbors() {
    neighbourhood neighbourhood;
    neighbourhood.determine_neighbors(node_infos);
}

/**
 * @fn void nodeGenerator::determine_neighbors()
 * @brief does what it says, no real checking and correcting thou !!
 * @return returns false if a file does not exist
 */
bool nodeGenerator::read_data_from_file() {
    // redo the struct
    if(redo) {return false;}
    // functional body
    bool return_value = false;
    std::filesystem::path file{file_name};
    if(std::filesystem::exists(file)) {
        std::ifstream data_file;
        std::string line;
        data_file.open(file_name);
        while(getline(data_file,line)) {
            // chop up the line to extract the information
            std::string delimiter = "|";
            size_t pos = 0;
            readBack_t chop = HANDLE;
            std::string token;
            // create a node
            auto node_point = new nodePoint_t;
            while ((pos = line.find(delimiter)) != std::string::npos) {
                token = line.substr(0, pos);
                read_back_switch_case(node_point,token,&chop);
                line.erase(0, pos + delimiter.length());
            }
            // last part is not choped off and still in the line string
            read_back_switch_case(node_point,line,&chop);
            node_infos.push_back(node_point);
        }
        return_value = true;
    }
    return return_value;
}
/**
 * @fn void nodeGenerator::read_back_switch_case(nodePoint_t* n, std::string& s, readBack_t* chop)
 * @brief function to translate data on files to data for the class
 * @param n
 * @param s
 * @param chop
 */
void nodeGenerator::read_back_switch_case(nodePoint_t* n, std::string& s, readBack_t* chop) {
    switch(*chop) {
    case HANDLE:
        n->handle = std::stoull(s); // hopefully correct
        *chop = TYPE;
        break;
    case TYPE:
        n->type = static_cast<nodeIdentifier_t>(std::stoi(s));
        *chop = BOUNDARY;
        break;
    case BOUNDARY:
        n->boundary = static_cast<boundaryType_t>(std::stoi(s));
        *chop = POSITION;
        break;
    case POSITION: {
        point_t point;
        n->position = point; // trick to set the array to the size of a point
        std::string delimiter = ", ";
        std::vector<std::string> split = split_string(s,delimiter);
        int i = 0;
        for(const auto& pos : split) {
            n->position(i) = std::stod(pos);
            i++;
        }
        *chop= LINKS;
        break;
    }
    case LINKS: {
        std::string delimiter = ";";
        std::string komma = ",";
        std::vector<std::string> split = split_string(s,delimiter);
        for(auto part : split) {
            std::vector<std::string> sub_split = split_string(part,komma);
            toLinks_t link;
            link.channel = std::stoi(sub_split.at(0));
            link.handle = std::stoull(sub_split.at(1));
            n->links.push_back(link);
        }
        // set the next state
        *chop = ERROR;
        break;
    }
    default:
        throw std::invalid_argument("Error while parsing");
    }
}
/**
 * @fn void nodeGenerator::write_data_to_file(bool write)
 * @brief does what it says check if params are set correctlly
 * @param write
 */
void nodeGenerator::write_data_to_file(bool write) {
    /// only write when given the task
    if(!write) {
        return;
    }
    std::ofstream  out;
    out.open(file_name);
    for(const auto i : node_infos) {
        // todo investigate if there is a better way to hand out data than in a file will have to do for now thou
        // general type info
        out << i->handle << " | "
            << i->type << " | "
            << i->boundary << " | " ;
        // position
        auto iter = i->position.begin();
        // unlike std containers no trailing end is added
        do {
            out << iter.operator*();
            iter++;
            if(iter != i->position.end()) {
                out <<  ", ";
            }
        }while(iter != i->position.end());
        out << " | ";
        // links + neighbours handles
        auto it = i->links.begin();
        // unlike std containers no trailing end is added
        while (it != i->links.end()) {
            out << it.operator*().channel << ", "
                << it.operator*().handle;
            if(it < i->links.end()-1) {
                out << " ; ";
            }
            it++;
        }
        out << std::endl;
    }
    out.close();
}

/**
 * @fn void nodeGenerator::board_creation(unsigned int size)
 * @brief create a board based on sizes given
 */
void nodeGenerator::board_creation(unsigned int size) {
    // make the interleaved positions
    std::vector<uint64_t> interleaved_positions;
    for(uint32_t i = 0; i < size; ++i) {
        for(uint32_t j=0; j < size; ++j) {
            uint64_t interleaved_position = bit_interleaving_2d(i,j);
            interleaved_positions.push_back(interleaved_position);
        }
    }
    // sort
    std::sort(interleaved_positions.begin(),interleaved_positions.end());
    // create a drawing board of nodes with side lengths size
    handle_t handle_counter = 1;
    for(auto p: interleaved_positions) {
        uint32_t x = bit_extraleaving_2d_x(p);
        uint32_t y = bit_extraleaving_2d_y(p);
        point_t point = {x,y};
        auto n = new nodePoint_t;
        n->handle = handle_counter;
        n->position = point;
        n->type = UNKNOWN;
        n->boundary = NO_BOUNDARY;
        // dont forget to increase the handle counter each time
        handle_counter++;
        node_infos.push_back(n);
    }
    /*
    handle_t handle_counter = 1;
    for(int i = 0; i < size; ++i) {
        for(int j = 0; j < size; ++j) {
            point_t point = {i,j};
            auto n = new nodePoint_t;
            n->handle = handle_counter;
            n->position = point;
            n->type = UNKNOWN;
            n->boundary = NO_BOUNDARY;
            // dont forget to increase the handle counter each time
            handle_counter++;
            node_infos.push_back(n);
        }
    }
     */
}

/**
 * @fn void nodeGenerator::check_nodes(handle_t* current)
 * @brief check if the node is inside the canvas defined by the straights and deletes the rest
 * @param current
 */
void nodeGenerator::check_nodes(handle_t* current) {
    straightGenerator straight(points);
    straight.init();
    std::vector<nodePoint_t*> reformed_nodes;
    for(auto n : node_infos) {
        if(!straight.node_inside(n)) {
            n->handle = *current;
            n->boundary = NO_BOUNDARY;
            n->type = WET;
            reformed_nodes.push_back(n);
            (*current)++;
        }
        else {
            // also delete element if not inside
            delete n;
        }
    }
    node_infos = reformed_nodes;
}

/**
 * @fn void nodeGenerator::add_boundary_nodes(handle_t* current)
 * @brief does what it says after sorting out everything we still have to add the boundary nodes
 * @param current
 */
void nodeGenerator::add_boundary_nodes(handle_t* current) {
    // lastlly add the boundary points
    for(auto p : points->boundary_points) {
        // std::cout << p << std::endl;
        // set all the variables except the links
        auto n = new nodePoint_t;
        n->handle = *current;
        n->position = p->point;
        n->type = DRY;
        n->boundary = p->type;
        // dont forget to increase the handle counter each time
        (*current)++;
        node_infos.push_back(n);
    }
}

/**
 * @fn
 * @brief
 */
void nodeGenerator::reduce_boundary_neighborhood() {
    int boundary_start = 0;
    for(auto n : node_infos) {
        if(n->type == DRY) {
            for(auto link : n->links) {
                handle_t partner_handle = link.handle;
                int link_channel = link.channel; // channel where the info is
                int from_channel = switch_link_dimensions(link_channel); // channel where it has to go
                long array_position = long(partner_handle) - 1;
                int channel_position = link_channel-1;
                // switchero
                node_infos.at(array_position)->links.at(channel_position).channel = from_channel;
                node_infos.at(array_position)->links.at(channel_position).handle = partner_handle;
            }
        }
        else {
            assert(n->links.size() == 8);
            boundary_start++;
        }
    }
    // delete because now unnecessary
    node_infos.erase(node_infos.begin() + boundary_start, node_infos.end());
}
/// public
/**
 * @fn void nodeGenerator::set_discovery_vector(vector_t set)
 * @brief set the 2D discovery vector, the function linear generation will use that vector during node discovery
 * @param set the discovery vector
 */
void nodeGenerator::set_discovery_vector(vector_t set) {
    discovery_vector = set;
}
/**
 * @fn void nodeGenerator::set_redo_save(bool r, bool s)
 * @brief
 * @param r redo
 * @param s safe
 */
void nodeGenerator::set_redo_save(bool r, bool s) {
    redo = r;
    save = s;
}

/**
 * @fn void nodeGenerator::init()
 * @brief initializes the node generator, if there are nodes given in the form of a stored_nodes_file, will use that
 * old legacy method
 */
void nodeGenerator::init() {
    if(!read_data_from_file()) {
        linear_generation();
        determine_neighbors();
        write_data_to_file(save);
    }
}

/**
 * @fn void nodeGenerator::init(unsigned int size)
 * @param size
 */
void nodeGenerator::init(unsigned int size) {
    if(!read_data_from_file()) {
        board_creation(size);
        handle_t handle_counter = 1;
        check_nodes(&handle_counter);
        add_boundary_nodes(&handle_counter);
        determine_neighbors();
        write_data_to_file(save);
    }
}

void nodeGenerator::init_fused(unsigned int size) {
    if(!read_data_from_file()) {
        board_creation(size);
        handle_t handle_counter = 1;
        check_nodes(&handle_counter);
        add_boundary_nodes(&handle_counter);
        determine_neighbors();
        reduce_boundary_neighborhood();
        write_data_to_file(save);
    }
}

/**
 * @fn void nodeGenerator::delete_node_infos()
 * @brief deletes the content of the node_infos vector
 */
void nodeGenerator::delete_node_infos() {
    // delete node info
    for(auto n : node_infos) {
        delete n;
    }
}