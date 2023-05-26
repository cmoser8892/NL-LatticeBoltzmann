# Defining the executable target
add_executable(06a main.cpp)

# Linking against everything
target_link_libraries(06a PUBLIC my_lib)
