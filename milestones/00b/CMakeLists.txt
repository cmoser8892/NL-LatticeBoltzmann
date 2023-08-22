# Defining the executable target
add_executable(00b main.cpp)

# Linking against everything
target_link_libraries(00b PUBLIC source_lib)
