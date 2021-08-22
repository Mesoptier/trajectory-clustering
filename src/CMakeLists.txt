file(GLOB_RECURSE HEADER_FILES *.h)

set(Eigen3_DIR "/usr/include/eigen-3.4-rc1/build")
find_package(Eigen3 3.3 REQUIRED NO_MODULE)

set(SOURCE_FILES
    geom.cpp
    io.cpp
    Curve.cpp
    IntegralFrechet/IntegralFrechet.cpp
    IntegralFrechet/Cell.cpp
    IntegralFrechet/metrics/L1.cpp
    IntegralFrechet/metrics/LInfinity_NoShortcuts.cpp
    SymmetricMatrix.cpp
    clustering/pam.cpp
    IntegralFrechet/MatchingBand.cpp
    simplification/imaiiri.cpp
    DTW/dtw.cpp
    cdtw/Polynomial.cpp
    cdtw/Interval.cpp
)

add_library(TrajectoryClustering_lib STATIC ${SOURCE_FILES} ${HEADER_FILES})
target_include_directories(TrajectoryClustering_lib PRIVATE ../src)
target_link_libraries(TrajectoryClustering_lib Eigen3::Eigen)