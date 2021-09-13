aux_source_directory(. SRC)
add_executable(hyrax_run ${SRC})
target_link_libraries(hyrax_run hyrax_lib mcl mclbn384_256)
add_test(NAME hyrax_run COMMAND hyrax_run)