cmake_minimum_required(VERSION 3.10)
project(code)

set(CMAKE_CXX_STANDARD 14)

find_package(Armadillo REQUIRED)
find_package(Boost REQUIRED)

set(SOURCE_FILES src/main.cpp src/FastMarchIntegralFrechet.cpp src/FastMarchIntegralFrechet.h src/Curve.h src/IntegralFrechet/Cell.cpp src/IntegralFrechet/Cell.h src/IntegralFrechet/Solver.cpp src/IntegralFrechet/Solver.h src/Vertex.h src/Edge.h)
add_executable(code ${SOURCE_FILES})

target_include_directories(
        code PRIVATE
        ${ARMADILLO_INCLUDE_DIRS}
        ${BOOST_INCLUDE_DIRS}
)
target_link_libraries(
        code
        ${ARMADILLO_LIBRARIES}
        ${BOOST_LIBRARIES}
)