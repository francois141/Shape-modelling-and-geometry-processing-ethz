cmake_minimum_required(VERSION 3.1)
include(FetchContent)
project(assignment6)

FetchContent_Declare(
        libigl
        GIT_REPOSITORY https://github.com/libigl/libigl.git
        GIT_TAG v2.4.0
)
FetchContent_MakeAvailable(libigl)
igl_include(glfw)
igl_include(imgui)

# Add your project files
file(GLOB SRC_FILES *.cpp)
add_executable(${PROJECT_NAME} ${SRC_FILES} util/Quaternion.h)
target_link_libraries(${PROJECT_NAME} PUBLIC igl::core igl::imgui igl::glfw)