cmake_minimum_required(VERSION 3.16)

project(ssm)

set(CMAKE_CXX_STANDARD 20)

add_library(ssm INTERFACE)
target_include_directories(ssm INTERFACE include)

include(FetchContent)

option(SSM_INCLUDE_CGAL "Include CGAL" OFF)
option(SSM_BUILD_APP "Build the sample GUI" OFF)
option(SSM_BUILD_TESTS "Build tests when building different projects" OFF)
option(SSM_BUILD_EXAMPLES "Build examples" ON)
option(SSM_BUILD_DOCS "Build documentation" OFF)

set(CGAL_FIND_COMPONENTS Core)

if(SSM_BUILD_EXAMPLES)
  list(APPEND CGAL_FIND_COMPONENTS Qt6)
endif()

if(SSM_INCLUDE_CGAL)
  FetchContent_Declare(
    cgal
    GIT_REPOSITORY https://github.com/CGAL/cgal.git
    GIT_TAG v6.0
  )

  FetchContent_MakeAvailable(cgal)

  include("${cgal_SOURCE_DIR}/CGALConfig.cmake")
else()
  find_package(CGAL 6.0 REQUIRED COMPONENTS ${CGAL_FIND_COMPONENTS})
endif()

FetchContent_Declare(
  thread_pool
  GIT_REPOSITORY https://github.com/bshoshany/thread-pool.git
  GIT_TAG v4.1.0
)

FetchContent_MakeAvailable(thread_pool)
target_include_directories(ssm INTERFACE ${thread_pool_SOURCE_DIR}/include)

target_link_libraries(ssm INTERFACE CGAL::CGAL)

if(SSM_BUILD_EXAMPLES)
  add_subdirectory(examples)
endif()

if((CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME OR SSM_BUILD_TESTS)
  AND BUILD_TESTING)
  enable_testing()
  add_subdirectory(tests)
endif()

if(SSM_BUILD_DOCS)
  find_package(Doxygen REQUIRED)

  doxygen_add_docs(
    ssm_docs
    include
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    COMMENT "Generating documentation with Doxygen"
  )
endif()
