aux_source_directory(. LIB_SRC)
list(FILTER LIB_SRC EXCLUDE REGEX "main*")

add_library(hyrax_lib ${LIB_SRC})
add_executable(hyrax_time main.cpp)
target_link_libraries(hyrax_time hyrax_lib mcl mclbn384_256)