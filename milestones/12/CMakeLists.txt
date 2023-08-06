# Defining the executable target
add_executable(12 main.cpp)

# Linking against everything
target_link_libraries(12 PUBLIC simulation_lib)
