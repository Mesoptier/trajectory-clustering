cmake_minimum_required(VERSION 3.10)
project(code)

set(CMAKE_CXX_STANDARD 17)

set(SOURCE_FILES main.cpp)
add_executable(ProjectName_run ${SOURCE_FILES})

# TODO: Rename directory to ProjectName_lib?
include_directories(src)
add_subdirectory(src)

target_link_libraries(ProjectName_run ProjectName_lib)

add_subdirectory(Google_tests)
