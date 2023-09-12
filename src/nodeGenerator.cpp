#include "nodeGenerator.h"
#include "helper_functions.h"
#include "neighborhood.h"
#include "functions.h"
#include <iostream>
#include <fstream>
#include <filesystem>
#include <sstream>

/**
 * Constructor, sets points.
 * @param p
 */
nodeGenerator::nodeGenerator(boundaryPointConstructor *p) {
    points = p;
    if(p != nullptr) {
        straight_surfaces = new straightGenerator(points);
        straight_surfaces->init();
    }
}

/**
 * Constructor, sets a surface.
 * @param s
 */
nodeGenerator::nodeGenerator(straightGenerator *s) {
    straight_surfaces = s;
    create_straight = false;
}

/**
 * Deletes all the info about the nodes.
 */
nodeGenerator::~nodeGenerator() {
    delete_node_infos();
    if(create_straight)
        delete straight_surfaces;
    //
    delete markers;
}

/**
 * Generates nodes until it hits a boundary useful to make simple 1d tests.
 */
void nodeGenerator::linear_generation() {
    handle_t handle_counter = 1;
    // go through the boundary points starting at a b point and go throu while still discovering new ones
    for(auto bs : points->boundary_structures) {
        for(auto p : bs->boundary_points) {
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
    }
    // finally add boundary nodes
    add_boundary_nodes(&handle_counter);
}

/**
 * Checks if the next point would be a boundary or outside.
 * @param p
 * @param check_point
 * @return true or false based if boundary/not also checks if still inside bounds
 */
bool nodeGenerator::check_other_boundary_hit(boundaryPoint_t* p,point_t &check_point) {
    // returns true if hit or outside false if not
    bool return_value = false;
    point_t check = check_point.base();
    // right now just does a linear search to check if point hit or not
    for(auto bs : points->boundary_structures) {
        for(auto point : bs->boundary_points) {
            if(compare_two_points(&check,&point->point)) {
                return_value = true;
                break;
            }
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
 * Hull for the neighbourhood class with its determine_neighbors call.
 */
void nodeGenerator::determine_neighbors() {
    neighbourhood neighbourhood;
    neighbourhood.determine_neighbors(node_infos);
}

/**
 * Does what it says, no real checking and correcting thou.
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
 * Function to translate data on files to data for the class.
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
 * Does what it says check if params are set correctly.
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
 * Create a board based on sizes given.
 */
void nodeGenerator::board_creation(unsigned int size) {
    size_canvas = size;
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
        n->boundary = INIT_NONE;
        // dont forget to increase the handle counter each time
        handle_counter++;
        // container updates
        node_infos.push_back(n);
        // add/setup the to be removed control
        to_be_removed.push_back(true);
    }
    fill_search();
}

/**
 * Check if the node is inside the canvas defined by the straights and sets those to not be deleted.
 */
void nodeGenerator::check_nodes_inside() {
    for(int i = 0; i < node_infos.size(); ++i) {
        auto n = node_infos[i];
        bool not_outside = straight_surfaces->node_inside_simple(n);
        if(!not_outside) {
            n->type = WET;
            n->boundary = NO_BOUNDARY;
            // should not be removed after all
            to_be_removed[i] = false;
        }
    }
}

/**
 * Checks the ibm relationship and flags the nodes.
 * @param range
 */
void nodeGenerator::check_nodes_ibm(double range) {
    // look through the markers and
    for(auto mark : markers->marker_points) {
        std::vector<handle_t> affected = rpkh.ranging_key_translation(*mark,range);
        for(auto handle : affected) {
            handle = handle - 1; // handle to the guy
            auto current_node = node_infos[handle];
            // check if inside
            if(current_node->boundary == NO_BOUNDARY) {
                current_node->type = WET;
                current_node->boundary = IBM_INNER;
            }
            // dont relabel
            else if(current_node->boundary == INIT_NONE) {
                current_node->type = WET;
                current_node->boundary = IBM_OUTER;
            }
            to_be_removed[handle] = false;
        }
    }
}

/**
 * Fills the range point key hash NL search.
 */
void nodeGenerator::fill_search() {
    handle_t current = 1;
    for(auto node : node_infos) {
        rpkh.fill_key(current,node->position);
        ++current;
    }
}

/**
 * Sets up the nodes for periodic boundaries.
 * @param t
 * @param a
 */
void nodeGenerator::check_nodes_periodic(kernelType_t t, long* a) {
    // need ibm range, to remove
    // creat markers for the surfaces affected
    for(int i = 0; i < 2; ++i) {
        // create and so
        periodic_marker[i] = new markerPoints();
        periodic_marker[i]->distribute_markers_periodic(straight_surfaces->surfaces[a[i]],IBM,t);
    }
    // first pass remove the unwanted node
    double range = kernel_id_to_lattice_search(KERNEL_C);
    // determine fluid direction from the middle nodes
    int run_variable = 0; // r
    for(auto pm : periodic_marker) {
        // find the middle marker
        auto s = double(pm->marker_points.size());
        auto middle = size_t(std::round(s/2.0));
        point_t midpoint = *pm->marker_points[middle];
        // look around
        std::vector<handle_t> found = rpkh.ranging_key_translation(midpoint,2.5);
        // addup everything that we found
        vector_t add_up = {0,0};
        for(auto handle: found) {
            handle -= 1;
            auto current_node = node_infos[handle];
            add_up += point_t(current_node->position) - midpoint;
        }
        // set the direction of the reference
        periodic_reference[run_variable] = -1 * add_up.normalized();
        // recenter
        periodic_reference[run_variable] = vector_to_cardinal(periodic_reference[run_variable]);
        // increment
        ++run_variable;
    }
    // remove excessive ibm nodes
    // loop over both
    run_variable = 0;
    for(auto pm : periodic_marker) {
        // loop over marker
        for(auto m : pm->marker_points) {
            // find affected nodes
            std::vector<handle_t> affected = rpkh.ranging_key_translation(*m,range);
            for(auto handle : affected) {
                // access to var
                handle -= 1;
                auto current_node = node_infos[handle];
                // make a vector to check
                vector_t checker = point_t(current_node->position) - *m;
                // reference goes out so if the vector from the position out +-90degrees its out
                if(check_plus_minus_90(&checker,&periodic_reference[run_variable])) {
                    // flag to be deleted
                    to_be_removed[handle] = true;
                }
            }
        }
        // increment var
        ++run_variable;
    }
    // reintroduce the needed periodic nodes ( i know this is extra work but this is a computer)
    // clarity is more important
    run_variable = 0;
    for(auto pm: periodic_marker) {
        for(auto m : pm->marker_points) {
            std::vector<handle_t> affected = rpkh.ranging_key_translation(*m,0.9);
            for(auto handle : affected) {
                handle -= 1;
                auto current_node = node_infos[handle];
                point_t node_point = current_node->position;
                if(compare_two_points(&node_point,m)) {
                    current_node->type = PERIODIC_CONNECT;
                    to_be_removed[handle] = false;
                }
            }
        }
    }
}

/**
 * Removes all the nodes that were not wanted in an abortion like process.
 * @param current
 */
void nodeGenerator::remove_unwanted_nodes(handle_t *current) {
    std::vector<nodePoint_t*> reformed_nodes;
    for(int i = 0; i < to_be_removed.size(); ++i) {
        auto control = to_be_removed[i];
        auto node = node_infos[i];
        // remove
        if(control) {
            delete node;
        }
        // keep
        else {
            // add to the reformed nodes
            node->handle = *current;
            // node->type = WET;
            reformed_nodes.push_back(node);
            // increment
            ++(*current);
        }
    }
    // overwrite
    node_infos = reformed_nodes;
}

/**
 * Does what it says after sorting out everything we still have to add the boundary nodes.
 * @param current
 */
void nodeGenerator::add_boundary_nodes(handle_t* current) {
    // lastlly add the boundary points
    for(auto bs : points->boundary_structures) {
        for(auto p : bs->boundary_points) {
            // std::cout << p << std::endl;
            // set all the variables except the links
            auto n = new nodePoint_t;
            n->handle = *current;
            n->position = p->point;
            n->type = p->dw;
            n->boundary = p->type;
            // dont forget to increase the handle counter each time
            (*current)++;
            node_infos.push_back(n);
        }
    }
}

/**
 * Function to delete the boundary nodes and write it directly in the bounce-back.
 */
void nodeGenerator::reduce_boundary_neighborhood() {
    // otherwise it just crashes
    int boundary_start = 0;
    for(auto n : node_infos) {
        if(n->type == DRY) {
            for(auto link : n->links) {
                handle_t partner_handle = link.handle;
                int link_channel = link.channel; // channel where the info is
                int from_channel = switch_link_dimensions(link_channel); // channel where it has to go
                long array_position = long(partner_handle) - 1;
                int channel_position = from_channel-1;
                // fall back assert
                assert(node_infos.size()>= array_position);
                assert(node_infos.at(array_position)->links.size() == (CHANNELS-1));
                // switch
                node_infos.at(array_position)->links.at(channel_position).channel = link_channel;
                node_infos.at(array_position)->links.at(channel_position).handle = partner_handle;
                check_and_set_reduced_neighborhood(array_position,n->boundary);
            }
        }
        else if(n->type == WET){
            //std::cerr << n->links.size() << std::endl;
            //std::cerr << n->position << std::endl;
            // assert(n->links.size() == 8);
            boundary_start++;
        }
        else {
            std::cerr << "node-generator: unknown node type" << std::endl;
        }
    }
    // delete because now unnecessary
    auto iter = node_infos.begin() + boundary_start;
    while(iter != node_infos.end()) {
        delete iter.operator*();
        ++iter;
    }
    // erase dangeling pointers
    node_infos.erase(node_infos.begin() + boundary_start, node_infos.end());
}

/**
 * Bumps up the the boundary condition of a wet node to represent the boundary.
 * @param array_position
 * @param b
 */
void nodeGenerator::check_and_set_reduced_neighborhood(handle_t array_position, boundaryType_t b) {
    // bumps up the boundary condition
    auto n = node_infos.at(array_position);
    auto boundary = n->boundary;
    // might disagree but the fall throughs of switch cases make this a bit more readable
    // if either bb or moving bb bump up
    // the others get their special treatment too
    switch(b) {
    case BOUNCE_BACK: {
        // fall through BB
        switch(boundary) {
        case NO_BOUNDARY:
            n->boundary = b;
            break;
        default :
            break;
        }
        break;
    }
    case BOUNCE_BACK_MOVING: {
        // fall through BB Moving
        switch(boundary) {
        case NO_BOUNDARY:
        case BOUNCE_BACK:
            n->boundary = b;
            break;
        default:
            break;
        }
        break;
    }
    case OPEN_INLET: {
        // fall through open inlet
        switch(boundary) {
        case NO_BOUNDARY:
        case BOUNCE_BACK:
        case BOUNCE_BACK_MOVING:
            n->boundary = b;
            break;
        default:
            break;
        }
        break;
    }
    default:
        break;
    }
}

/**
 * Approximates the outer boundary with a bb boundary.
 * @note makes sure every nodes has 8 links
 */
void nodeGenerator::fill_neighborhood_holes() {
    int counter = 0;
    for(auto ni : node_infos) {
        // find the boarder nodes
        int link_channel = 1;
        if(ni->links.size() < 8) {
            counter++;
            for(int pos = 0; pos < 8; ++pos) {
                toLinks_t link;
                if(pos < ni->links.size()) {
                    link = ni->links[pos];
                }
                else {
                    link.channel = 0;
                    link.handle = 0;
                }
                // expected channel
                if(link.channel == link_channel) {
                    // nop
                }
                else {
                    // create a new link and insert it
                    toLinks_t new_link;
                    new_link.handle = ni->handle;
                    int linked_channel = switch_link_dimensions(link_channel);
                    new_link.channel = linked_channel;
                    ni->links.insert(ni->links.begin()+pos ,new_link);
                }
                ++link_channel;
            }
        }
    }
}


void nodeGenerator::connect_periodic_boundary() {
    // it is important to keep this distinction! to be able to resolve it with the reference!
    periodicBundles bundle(periodic_marker[0],periodic_marker[1]);
    bundle.connect();
    // rebuild rpkh
    rpkh.clear();
    fill_search();
    for( auto ass : bundle.associations) {
        // resolve the association -> find both nodes
        handle_t inlet_location = ass.first;
        handle_t outlet_location = ass.second;
        // makers
        point_t markers_relevant[2] = {*periodic_marker[0]->marker_points[inlet_location],
                                       *periodic_marker[1]->marker_points[outlet_location]};
        // find the nodes
        std::vector<handle_t> affected_inlet
            = rpkh.ranging_key_translation(markers_relevant[0],0.9);
        std::vector<handle_t> affected_outlet
            = rpkh.ranging_key_translation(markers_relevant[1],0.9);
        if((affected_inlet.size() != 1) || (affected_outlet.size() != 1)){
            std::cerr << "Bad inlet/outlet" << std::endl;
        }
        // take the first two
        nodePoint_t  * inlet_outlet_nodes[2] = {node_infos[affected_inlet[0]-1],
                                                node_infos[affected_outlet[0]-1]};
        // for those two resolve the neighborhood list !
        // todo this feels and looks like italian pasta
        for(int i = 0; i < 2; ++i) {
            // tested with the little_things suit
            int other = (i+1)%2;
            // init parameters
            auto self = inlet_outlet_nodes[i];
            auto partners_point = markers_relevant[other];
            vector_t set_vector = periodic_reference[i];
            // swap x and y
            vector_t manipulator = {set_vector.y(),set_vector.x()};
            for (int j = 0; j < 3; ++j) {
                // manipulate partner and set vector
                // we need a combined manipulator 0 -> 0, 1 -> 1, 2 -> -1
                // i know this is dumb...
                int first_manipulator = 0;
                if(j > 0) {
                    first_manipulator = 1;
                }
                int second_manipulator = 1;
                if(j == 2) {
                    second_manipulator = -1;
                }
                // calculate the velocity set vector
                vector_t setter_vector = set_vector + first_manipulator*second_manipulator*manipulator;
                // calculate the position of the partner
                point_t partners_point_position = partners_point + first_manipulator*second_manipulator*manipulator;
                // set the relevant parameters in the link fields
                set_periodic_boundary(self,&partners_point_position,&setter_vector);
            }
        }
    }
}

void nodeGenerator::set_periodic_boundary(nodePoint_t* self,
                                          point_t* partner_position,
                                          vector_t* set) {
    // find the handle of the partner
    std::vector<handle_t> search = rpkh.ranging_key_translation(*partner_position,0.9);
    int channel_id = index_of_velocity_set(*set);
    // if valid 
    if(search.size() == 1) {
        handle_t valid_handle = search[0];
        self->links[channel_id-1].channel = channel_id;
        self->links[channel_id-1].handle = valid_handle;
    }
}

/// public
/**
 * Set the 2D discovery vector, the function linear generation will use that vector during node discovery.
 * @param set the discovery vector
 */
void nodeGenerator::set_discovery_vector(vector_t set) {
    discovery_vector = set;
}

/**
 * Sets the redo save, weather or not we save to text or not.
 * @param r redo
 * @param s safe
 */
void nodeGenerator::set_redo_save(bool r, bool s) {
    redo = r;
    save = s;
}

/**
 * Initializes the node generator, if there are nodes given in the form of a stored_nodes_file, will use that.
 * old legacy method
 * @note used for liniar structures
 */
void nodeGenerator::init() {
    if(!read_data_from_file()) {
        linear_generation();
        determine_neighbors();
        write_data_to_file(save);
    }
}

/**
 * Init.
 * @note Inits a bounce back structure and leaves the boundary bb nodes in.
 * @param size
 */
void nodeGenerator::init(unsigned int size) {
    if(!read_data_from_file()) {
        board_creation(size);
        handle_t handle_counter = 1;
        check_nodes_inside();
        remove_unwanted_nodes(&handle_counter);
        add_boundary_nodes(&handle_counter);
        determine_neighbors();
        write_data_to_file(save);
    }
}

/**
 * Fused init, also reduces total nodes by removing boundaries.
 * @note Same as init(), but the bounce back nodes are optimized away.
 * @param size canvas size
 */
void nodeGenerator::init_fused(unsigned int size) {
    if(!read_data_from_file()) {
        board_creation(size);
        handle_t handle_counter = 1;
        check_nodes_inside();
        remove_unwanted_nodes(&handle_counter);
        add_boundary_nodes(&handle_counter);
        determine_neighbors();
        reduce_boundary_neighborhood();
        write_data_to_file(save);
    }
}

/**
 * Inits the simulation based on a given surface.
 * @param size
 * @param range
 */
void nodeGenerator::init_surface(unsigned int size, double range) {
    // correct the range to
    if(!read_data_from_file()) {
        // create the drawing canvas
        board_creation(size);
        handle_t handle_counter = 1;
        // generate the markers
        if(markers == nullptr) {
            markers = new markerPoints(straight_surfaces);
            markers->distribute_markers();
        }
        // check and reform nodes
        check_nodes_inside();
        check_nodes_ibm(range);
        remove_unwanted_nodes(&handle_counter);
        determine_neighbors();
        write_data_to_file(save);
    }
}

/**
 * Tests out weather or not adding an additional bb outside improves stability (it kinda doesnt)
 * @param size
 * @param range
 * @param marker_range
 */
void nodeGenerator::init_surface_return(unsigned int size, kernelType_t type,  double marker_range) {
    double range = kernel_id_to_lattice_search(type);
    bool periodic_id = false;
    // correct the range to
    if(!read_data_from_file()) {
        // create the drawing canvas
        board_creation(size);
        handle_t handle_counter = 1;
        // generate the markers
        if(markers == nullptr) {
            markers = new markerPoints(straight_surfaces, marker_range);
            markers->distribute_markers();
        }
        // check and reform nodes
        check_nodes_inside();
        check_nodes_ibm(range);
        // periodic handeling
        // loop over the surface to find the periodic ones
        long holder[2] = {-1,-1};
        int i = 0;
        for(int s = 0; s < straight_surfaces->surfaces.size(); ++s) {
            auto surface = straight_surfaces->surfaces[s];
            if(surface->type == PERIODIC) {
                holder[i] = s;
                ++i;
                if(i >2) {
                    std::cerr << "To many periodic surfaces" << std::endl;
                }
            }
        }
        // check valid
        // todo periodics fkt do something wierd
        if((holder[0] >= 0) && (holder[2] >= 0)) {
            std::cout << "periodics" << std::endl;
            // check_nodes_periodic(type,holder);
            periodic_id = true;
        }
        // remove nodes
        remove_unwanted_nodes(&handle_counter);
        // neighborhood
        determine_neighbors();
        fill_neighborhood_holes();
        // modify for periodic
        if(periodic_id) {
            // create associations between inlet and outlet
            // rewrite the connections
            // connect_periodic_boundary();
        }
        write_data_to_file(save);
    }
}

/**
 * Deletes the content of the node_infos vector.
 */
void nodeGenerator::delete_node_infos() {
    // delete node info
    for(auto n : node_infos) {
        delete n;
    }
    // get rid of dangeling pointers
    node_infos.clear();
}

/**
 * Simple visualizer node points, if points where given uses that.
 * @param size
 */
void nodeGenerator::visualize_2D_nodes() {
    flowfield_t output;
    if(points != nullptr)
        output.setZero(std::floor(points->size.x()),std::floor(points->size.y()));
    else
        output.setZero(size_canvas,size_canvas);
    for(auto b : node_infos) {
        ++output(int(b->position.x()),int(b->position.y()));
    }
    std::cout << "Nodes allocated" << std::endl;
    std::cout << output << std::endl << std::endl;
}

/**
 * Visualize 2d nodes based on the label.
 * @param t
 */
void nodeGenerator::visualize_2D_nodes_labels(boundaryType_t t) {
    // setup the output field
    flowfield_t output;
    if(points != nullptr)
        output.setZero(std::floor(points->size.x()),std::floor(points->size.y()));
    else
        output.setZero(size_canvas,size_canvas);
    // give out based on th tag
    for(auto b : node_infos) {
        if(b->boundary == t)
            ++output(int(b->position.x()),int(b->position.y()));
    }
    std::cout << "Nodes with that tag: " << t << std::endl;
    std::cout << output << std::endl;
}

/**
 * Write out the nodes in file format useful for visualization.
 * @param t
 * @param write_file
 */
void nodeGenerator::write_out_nodes(boundaryType_t t, bool write_file) {
    // setup sizes
    flowfield_t output;
    if(points != nullptr)
        output.setZero(std::floor(points->size.x()),std::floor(points->size.y()));
    else
        output.setZero(size_canvas,size_canvas);
    // give ou the data
    for(auto b : node_infos) {
        if(b->boundary == t)
            ++output(int(b->position.x()),int(b->position.y()));
    }
    // condense together
    std::stringstream filename;
    filename<<  "node_type_file_" << t;
    //
    write_flowfield_data(&output, filename.str(),write_file);
}
