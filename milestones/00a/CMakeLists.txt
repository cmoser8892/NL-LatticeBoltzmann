# Defining the executable target
add_executable(00a main.cpp)

# Linking against everything
target_link_libraries(00a PUBLIC source_lib)
