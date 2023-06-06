# Defining the executable target
add_executable(09 main.cpp)

# Linking against everything
target_link_libraries(09 PUBLIC my_lib)
