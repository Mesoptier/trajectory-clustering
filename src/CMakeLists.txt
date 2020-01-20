file(GLOB_RECURSE HEADER_FILES *.h)

set(SOURCE_FILES
    geom.cpp
    io.cpp
    Curve.cpp
    IntegralFrechet/IntegralFrechet.cpp
    IntegralFrechet/Cell.cpp
    IntegralFrechet/metrics/L1.cpp
    IntegralFrechet/metrics/LInfinity_NoShortcuts.cpp
)

add_library(TrajectoryClustering_lib STATIC ${SOURCE_FILES} ${HEADER_FILES})
target_include_directories(TrajectoryClustering_lib PRIVATE ../src)