cmake_minimum_required(VERSION 3.10)
project(code)

set(CMAKE_CXX_STANDARD 14)

# TODO: clean-up these external libraries
include_directories(
        "C:/lib/armadillo-9.800.2/include"
)

link_libraries(
        "C:/lib/armadillo-9.800.2/examples/lib_win64/blas_win64_MT.lib"
        "C:/lib/armadillo-9.800.2/examples/lib_win64/lapack_win64_MT.lib"
)

set(SOURCE_FILES main.cpp)
add_executable(ProjectName_run ${SOURCE_FILES})

# TODO: Rename directory to ProjectName_lib?
include_directories(src)
add_subdirectory(src)

target_link_libraries(ProjectName_run ProjectName_lib)

add_subdirectory(Google_tests)
