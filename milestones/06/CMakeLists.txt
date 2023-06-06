# Defining the executable target
add_executable(06 main.cpp)

# Linking against everything
target_link_libraries(06 PUBLIC simulation_lib)
