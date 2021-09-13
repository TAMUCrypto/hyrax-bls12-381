cmake_minimum_required(VERSION 3.10)
project(hyrax_p224)
set(CMAKE_CXX_STANDARD 14)

include_directories(src)
include_directories(3rd/mcl/include)
include_directories(..)

add_subdirectory(src)
add_subdirectory(test)
add_subdirectory(3rd/mcl)
