file(GLOB_RECURSE HEADER_FILES *.h)

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
    clustering/pam_with_simplifications.cpp
    IntegralFrechet/MatchingBand.cpp
    clustering/clustering_algs.cpp
    curve_simplification.cpp
    clustering/center_algs.cpp
    simplification/imaiiri.cpp
    Frechet/frechet_light.cpp
    Frechet/frechet_matching.cpp
    filter.cpp
    clustering/center_clustering_algs.cpp
    DTW/dtw.cpp
    clustering/center_update.cpp
    random.cpp
)

add_library(TrajectoryClustering_lib STATIC ${SOURCE_FILES} ${HEADER_FILES})
target_include_directories(TrajectoryClustering_lib PRIVATE ../src)
target_link_libraries (TrajectoryClustering_lib Eigen3::Eigen)
