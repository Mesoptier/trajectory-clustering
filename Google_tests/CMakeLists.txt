project(Google_tests)
add_subdirectory(googletest)
include_directories(${gtest_SOURCE_DIR}/include ${gtest_SOURCE_DIR})

add_executable(Google_Tests_run geomTest.cpp CellTest.cpp)

target_link_libraries(Google_Tests_run TrajectoryClustering_lib)
target_link_libraries(Google_Tests_run gtest gtest_main)