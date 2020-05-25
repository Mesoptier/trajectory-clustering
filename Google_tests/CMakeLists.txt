project(Google_tests)
add_subdirectory(googletest)
include_directories(${gtest_SOURCE_DIR}/include ${gtest_SOURCE_DIR})

add_executable(Google_Tests_run
    geomTest.cpp
    ComputeMatchingTest.cpp
    NaiveLowerEnvelopeTest.cpp
    PolynomialTest.cpp
    CellTest.cpp
    CDTW_L1_L2Squared_L1.cpp
)

target_link_libraries(Google_Tests_run TrajectoryClustering_lib)
target_link_libraries(Google_Tests_run gtest gtest_main)