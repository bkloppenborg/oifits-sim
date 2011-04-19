cmake_minimum_required(VERSION 2.6) 
project(oifits-sim)

# Call CMake on the subdirectories

SET(EXECUTABLE_OUTPUT_PATH ${CMAKE_SOURCE_DIR}/bin)

add_subdirectory(lib)
add_subdirectory(src)
