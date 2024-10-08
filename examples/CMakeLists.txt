macro(package_add_example NAME)
    # create an executable in which the tests will be stored
    add_executable(${NAME} ${ARGN})

    # link the Google test infrastructure, mocking library, and a default main function to
    # the test executable.  Remove g_test_main if writing your own main function.
    target_link_libraries(${NAME} PRIVATE ssm CGAL::CGAL_Basic_viewer)
    target_include_directories(${NAME} PRIVATE ${PROJECT_SOURCE_DIR}/examples)
    set_target_properties(${NAME} PROPERTIES FOLDER examples)
endmacro()

add_subdirectory(SSM_restricted_voronoi_diagram)
