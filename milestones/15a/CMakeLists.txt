# Defining the executable target
add_executable(15a main.cpp)

# Linking against everything
target_link_libraries(15a PUBLIC simulation_lib)
