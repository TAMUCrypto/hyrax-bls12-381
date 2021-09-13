cmake_minimum_required(VERSION 3.10)
project(hyrax_p224)
set(CMAKE_CXX_STANDARD 14)

include_directories(./src)
include_directories(./lib)
include_directories(..)

add_subdirectory(src)
add_subdirectory(test)
