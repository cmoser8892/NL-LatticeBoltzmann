# Defining the executable target
add_executable(11b main.cpp)

# Linking against everything
target_link_libraries(11b PUBLIC simulation_lib)
