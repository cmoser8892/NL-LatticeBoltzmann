# Defining the executable target
add_executable(05 main.cpp)

# Linking against our MD code (propagates all dependencies)
target_link_libraries(05 PUBLIC simulation_lib)
