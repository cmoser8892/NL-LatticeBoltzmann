# Defining the executable target
add_executable(00 main.cpp)

# Linking against everything
target_link_libraries(00 PUBLIC simulation_lib)
