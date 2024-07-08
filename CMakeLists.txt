cmake_minimum_required(VERSION 3.15)

project(ssm)

set(CMAKE_CXX_STANDARD 20)

add_library(ssm INTERFACE)
target_include_directories(ssm INTERFACE include)

find_package(Eigen3 3.3 REQUIRED NO_MODULE)
find_package(CGAL REQUIRED)
target_link_libraries(ssm INTERFACE Eigen3::Eigen CGAL::CGAL)

option(SSM_BUILD_APP "Build the sample GUI" ON)
option(SSM_INCLUDE_AUTODIFF "Include autodiff" OFF)
option(SSM_ENABLE_AUTODIFF "Enable autodiff" OFF)

include(FetchContent)
if(SSM_BUILD_APP)
    add_executable(ssm_app src/main.cpp)

    FetchContent_Declare(
      libigl
        GIT_REPOSITORY https://github.com/libigl/libigl.git
        GIT_TAG        v2.5.0
    )

    # set(LIBIGL_WITH_CGAL ON)
    FetchContent_MakeAvailable(libigl)
    igl_include(glfw)
    igl_include(imgui)
    igl_include(stb)

    target_link_libraries(ssm_app PRIVATE ssm igl::core igl::opengl igl::glfw igl::imgui igl::stb)
endif()

if(SSM_INCLUDE_AUTODIFF)
    FetchContent_Declare(
      autodiff
        GIT_REPOSITORY https://github.com/autodiff/autodiff.git
        GIT_TAG        v1.1.2
    )
    
    set(AUTODIFF_BUILD_TESTS OFF)
    set(AUTODIFF_BUILD_EXAMPLES OFF)
    set(AUTODIFF_BUILD_PYTHON OFF)
    set(AUTODIFF_BUILD_DOCS OFF)

    FetchContent_MakeAvailable(autodiff)
    target_link_libraries(ssm INTERFACE autodiff::autodiff)
endif()

if(SSM_ENABLE_AUTODIFF)
    target_compile_definitions(ssm INTERFACE SSM_ENABLE_AUTODIFF)
endif()
