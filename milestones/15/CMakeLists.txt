# Defining the executable target
add_executable(15 main.cpp)

# Linking against everything
target_link_libraries(15 PUBLIC simulation_lib)
