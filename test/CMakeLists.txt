aux_source_directory(. SRC)
add_executable(hyrax_run ${SRC})
target_link_libraries(hyrax_run hyrax_lib gmp)
add_test(NAME hyrax_run COMMAND hyrax_run)