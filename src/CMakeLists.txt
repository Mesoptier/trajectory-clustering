file(GLOB_RECURSE HEADER_FILES *.h)

set(SOURCE_FILES
    Curve.cpp
    filter.cpp
    geom.cpp
    io.cpp
    SymmetricMatrix.cpp
    clustering/center_algs.cpp
    clustering/center_clustering_algs.cpp
    clustering/clustering_algs.cpp
    clustering/pam.cpp
    clustering/pam_with_simplifications.cpp
    DTW/dtw.cpp
    Frechet/frechet_light.cpp
    Frechet/frechet_matching.cpp
    IntegralFrechet/Cell.cpp
    IntegralFrechet/IntegralFrechet.cpp
    IntegralFrechet/MatchingBand.cpp
    IntegralFrechet/metrics/L1.cpp
    IntegralFrechet/metrics/LInfinity_NoShortcuts.cpp
    simplification/agarwal.cpp
    simplification/imaiiri.cpp
)

add_library(TrajectoryClustering_lib STATIC ${SOURCE_FILES} ${HEADER_FILES})
target_include_directories(TrajectoryClustering_lib PUBLIC ./)
