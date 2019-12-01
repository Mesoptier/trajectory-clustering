project(ProjectName)

set(HEADER_FILES
        Curve.h
        IntegralFrechet/Cell.h
        IntegralFrechet/Solver.h
        Vertex.h
        Edge.h
        util.h
        IntegralFrechet/metrics.h
        geom.h
        expressionml.h
        )

set(SOURCE_FILES
        IntegralFrechet/Cell.cpp
        IntegralFrechet/Solver.cpp
        geom.cpp
        )

add_library(ProjectName_lib STATIC ${SOURCE_FILES} ${HEADER_FILES})
