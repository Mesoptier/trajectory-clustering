cmake_minimum_required(VERSION 3.10)
project(code)

set(CMAKE_CXX_STANDARD 14)

find_package(Armadillo REQUIRED)
include_directories(${ARMADILLO_INCLUDE_DIRS})

set(SOURCE_FILES main.cpp)

add_executable(code ${SOURCE_FILES})
target_link_libraries(code ${ARMADILLO_LIBRARIES})