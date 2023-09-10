# Defining the executable target
add_executable(14 main.cpp)

# Linking against everything
target_link_libraries(14 PUBLIC simulation_lib)
