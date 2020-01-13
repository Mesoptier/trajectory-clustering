project(ProjectName)

set(HEADER_FILES
        Curve.h
        Edge.h
        util.h
        IntegralFrechet/metrics.h
        id.h
        geom.h
        expressionml.h
        io.h
        a_star.h
        IntegralFrechet/IntegralFrechet.h
        IntegralFrechet/Cell.h
        )

set(SOURCE_FILES
        geom.cpp
        io.cpp
        Curve.cpp
        IntegralFrechet/IntegralFrechet.cpp
        IntegralFrechet/Cell.cpp
        )

add_library(ProjectName_lib STATIC ${SOURCE_FILES} ${HEADER_FILES})
