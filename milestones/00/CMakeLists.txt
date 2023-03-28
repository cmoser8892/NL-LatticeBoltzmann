# Defining the executable target
add_executable(00 main.cpp)

# Linking against our MD code (propagates all dependencies)
target_link_libraries(00 PUBLIC my_lib)
