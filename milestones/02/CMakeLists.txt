# Defining the executable target
add_executable(02 main.cpp)

# Linking against our MD code (propagates all dependencies)
target_link_libraries(02 PUBLIC my_lib)