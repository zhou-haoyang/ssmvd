cmake_minimum_required(VERSION 3.5)

project(ssm)

set(CMAKE_CXX_STANDARD 20)

add_library(ssm INTERFACE)
target_include_directories(ssm INTERFACE include)

include(FetchContent)

option(SSM_INCLUDE_CGAL "Include CGAL" OFF)
option(SSM_BUILD_APP "Build the sample GUI" OFF)
option(SSM_BUILD_TESTS "Build tests when building different projects" OFF)
option(SSM_BUILD_EXAMPLES "Build examples" ON)

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

target_link_libraries(ssm INTERFACE CGAL::CGAL)

if(SSM_BUILD_APP)
  add_executable(ssm_app src/main.cpp)

  FetchContent_Declare(
    libigl
    GIT_REPOSITORY https://github.com/libigl/libigl.git
    GIT_TAG v2.5.0
  )

  # set(LIBIGL_WITH_CGAL ON)
  FetchContent_MakeAvailable(libigl)
  igl_include(glfw)
  igl_include(imgui)
  igl_include(stb)

  find_package(Eigen3 3.3 REQUIRED NO_MODULE)

  target_link_libraries(ssm_app PRIVATE ssm igl::core igl::opengl igl::glfw igl::imgui igl::stb Eigen3::Eigen)
endif()

if(SSM_BUILD_EXAMPLES)
  add_subdirectory(examples)
endif()

if((CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME OR SSM_BUILD_TESTS)
  AND BUILD_TESTING)
  enable_testing()
  add_subdirectory(tests)
endif()
