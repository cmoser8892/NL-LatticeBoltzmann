# Defining the executable target
add_executable(11a main.cpp)

# Linking against everything
target_link_libraries(11a PUBLIC simulation_lib)
