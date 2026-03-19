cmake_minimum_required(VERSION 3.10)
project(GraphPartitioner)

set(CMAKE_CXX_STANDARD 20)

include_directories(external/maxflow)
add_subdirectory(external/maxflow)

add_executable(partitioner main.cpp partitioner.cpp)
target_include_directories(partitioner PUBLIC datastructures external)
target_link_libraries(partitioner maxflow)
