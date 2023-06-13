# Defining the executable target
add_executable(07 main.cpp)

# Linking against everything
target_link_libraries(07 PUBLIC simulation_lib)
