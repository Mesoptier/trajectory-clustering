cmake_minimum_required(VERSION 3.12)
project(trajectory-clustering)

set(CMAKE_CXX_STANDARD 17)

add_subdirectory(src)
add_subdirectory(Google_tests)

# Main executable
add_executable(TrajectoryClustering_run main.cpp)

target_link_libraries(TrajectoryClustering_run PRIVATE TrajectoryClustering_lib)