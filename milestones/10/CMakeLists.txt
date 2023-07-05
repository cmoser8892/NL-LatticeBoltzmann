# Defining the executable target
add_executable(10 main.cpp)

# Linking against everything
target_link_libraries(10 PUBLIC simulation_lib)
