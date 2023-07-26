# Defining the executable target
add_executable(11 main.cpp)

# Linking against everything
target_link_libraries(11 PUBLIC simulation_lib)
