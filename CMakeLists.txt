cmake_minimum_required(VERSION 3.10)
project(code)

set(CMAKE_CXX_STANDARD 14)

find_package(Armadillo REQUIRED)
find_package(Boost REQUIRED)

set(SOURCE_FILES src/main.cpp src/FastMarchCDTW.cpp src/FastMarchCDTW.h src/Curve.h)
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