# List of headers
set(MY_HEADERS
        one_step_simulation.h
        two_step_simulation.h
        simulation.h)

# List of implementation files
set(MY_CPP
        one_step_simulation.cpp
        two_step_simulation.cpp
        simulation.cpp)

# Create a static library
add_library(simulation_lib STATIC ${MY_CPP} ${MY_HEADERS})

# Include the current directory as INTERFACE
# This will propagate the include to any target linking against the lib
target_include_directories(simulation_lib INTERFACE ${CMAKE_CURRENT_SOURCE_DIR})

# Link against Eigen (automatically sets up the includes)
target_link_libraries(simulation_lib PUBLIC source_lib)