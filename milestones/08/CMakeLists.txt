# Defining the executable target
add_executable(08 main.cpp)

# Linking against everything
target_link_libraries(08 PUBLIC simulation_lib)
