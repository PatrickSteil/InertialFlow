cmake_minimum_required(VERSION 3.10)
project(GraphPartitioner)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(THREADS_PREFER_PTHREAD_FLAG ON)
find_package(Threads REQUIRED)

include_directories(external/maxflow)
add_subdirectory(external/maxflow)

add_executable(partitioner main.cpp partitioner.cpp)
target_include_directories(partitioner PUBLIC datastructures external)
target_link_libraries(partitioner maxflow Threads::Threads)

option(INERTIALFLOW_BUILD_TESTS "Build the InertialFlow unit tests" ON)
if(INERTIALFLOW_BUILD_TESTS)
  enable_testing()

  add_executable(inertialflow_tests
    tests/test_main.cpp
    tests/test_generation_checker.cpp
    tests/test_partitioner.cpp
    partitioner.cpp
  )
  target_include_directories(inertialflow_tests PUBLIC external tests ${CMAKE_SOURCE_DIR})
  target_link_libraries(inertialflow_tests maxflow Threads::Threads)

  include(CTest)
  add_test(NAME inertialflow_tests COMMAND inertialflow_tests)
endif()
