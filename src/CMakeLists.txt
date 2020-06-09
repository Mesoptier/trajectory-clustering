file(GLOB_RECURSE HEADER_FILES *.h)

set(SOURCE_FILES
    classification_experiment.cpp
    Curve.cpp
    distance_functions.cpp
    geom.cpp
    simplification_experiment.cpp
    synthetic_curves.cpp
    clustering/center_algs.cpp
    clustering/center_clustering_algs.cpp
    clustering/clustering_algs.cpp
    clustering/pam.cpp
    DTW/dtw.cpp
    Frechet/filter.cpp
    Frechet/frechet_abstract.cpp
    Frechet/frechet_light.cpp
    Frechet/frechet_matching.cpp
    IntegralFrechet/Cell.cpp
    IntegralFrechet/IntegralFrechet.cpp
    IntegralFrechet/MatchingBand.cpp
    IntegralFrechet/metrics/L1.cpp
    IntegralFrechet/metrics/LInfinity_NoShortcuts.cpp
    simplification/agarwal.cpp
    simplification/imaiiri.cpp
    utils/CurveSimpMatrix.cpp
    utils/io.cpp
    utils/SymmetricMatrix.cpp
)

add_library(TrajectoryClustering_lib STATIC ${SOURCE_FILES} ${HEADER_FILES})
target_include_directories(TrajectoryClustering_lib PUBLIC ./)
