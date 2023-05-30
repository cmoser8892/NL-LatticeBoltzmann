# Defining the executable target
add_executable(06b main.cpp)

# Linking against everything
target_link_libraries(06b PUBLIC my_lib)
