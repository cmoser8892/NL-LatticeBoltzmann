# List of headers
set(MY_HEADERS
        simulation.h types.h functions.h helper_functions.h boundary_point_generator.h nodeGenerator.h straight.h neighborhood.h)

# List of implementation files
set(MY_CPP
        simulation.cpp functions.cpp helper_functions.cpp boundary_point_generator.cpp node.h node.cpp nodeGenerator.cpp straight.cpp neighborhood.cpp)

# Create a static library
add_library(my_lib STATIC ${MY_CPP} ${MY_HEADERS})

# Include the current directory as INTERFACE
# This will propagate the include to any target linking against the lib
target_include_directories(my_lib INTERFACE ${CMAKE_CURRENT_SOURCE_DIR})

# Link against Eigen (automatically sets up the includes)
target_link_libraries(my_lib PUBLIC Eigen3::Eigen)

# Set up MPI includes and library linking
# This also propagates to further targets
if (MPI_FOUND)
  target_compile_options(my_lib PUBLIC -Wall -Wextra -Wpedantic)
  target_compile_definitions(my_lib INTERFACE USE_MPI)
  target_include_directories(my_lib PUBLIC ${MPI_CXX_INCLUDE_DIRS} SYSTEM)
  target_link_libraries(my_lib PUBLIC MPI::MPI_CXX)
endif()
