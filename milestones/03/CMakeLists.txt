# Defining the executable target
add_executable(03 main.cpp)

# Linking against our MD code (propagates all dependencies)
target_link_libraries(03 PUBLIC my_lib)