include(FetchContent)
FetchContent_Declare(
    googletest
    URL https://github.com/google/googletest/archive/03597a01ee50ed33e9dfd640b249b4be3799d395.zip
)

# For Windows: Prevent overriding the parent project's compiler/linker settings
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)

include(GoogleTest)

add_executable(
    ssm_test
    SSM_voronoi_diagram/Simple_cartesian_find_intervals.cpp
    SSM_voronoi_diagram/Simple_cartesian_random.cpp
)

# link the Google test infrastructure, mocking library, and a default main function to
# the test executable.  Remove g_test_main if writing your own main function.
target_link_libraries(ssm_test ssm gtest gmock gtest_main)
target_include_directories(ssm_test PRIVATE ${PROJECT_SOURCE_DIR}/tests)
target_compile_definitions(ssm_test PRIVATE CGAL_DEBUG)

# gtest_discover_tests replaces gtest_add_tests,
# see https://cmake.org/cmake/help/v3.10/module/GoogleTest.html for more options to pass to it
gtest_discover_tests(ssm_test

    # set a working directory so your project root so that you can find test data via paths relative to the project root
    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
    PROPERTIES VS_DEBUGGER_WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}"
)
set_target_properties(ssm_test PROPERTIES FOLDER tests)
