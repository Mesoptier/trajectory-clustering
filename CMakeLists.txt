cmake_minimum_required(VERSION 3.13)
project(trajectory-clustering LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_EXTENSIONS OFF)

# CMake sets -g by default for Debug build
set(DEBUG_OPTIONS "-Og")
set(RELEASE_OPTIONS "-O3")

option(ASAN "Use address sanitizer" Off)
option(UBSAN "Use UB sanitizer" Off)
option(PROF "Set up for use with gprof" Off)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()
# Debug, Release, RelWithDebugInfo, MinSizeRel

add_compile_options("$<IF:$<CONFIG:DEBUG>,${DEBUG_OPTIONS},${RELEASE_OPTIONS}>")

if(PROF)
    add_compile_options("-pg")
    add_link_options("-pg")
endif()

if(NOT MSVC AND (ASAN OR UBSAN))
    set(SANITIZER "$<IF:$<BOOL:${ASAN}>,address,undefined>")
    set(SAN_OPTS "-fno-omit-frame-pointer" "-fsanitize=${SANITIZER}")
    add_compile_options("$<$<CONFIG:DEBUG>:${SAN_OPTS}>")
    add_link_options("$<$<CONFIG:DEBUG>:${SAN_OPTS}>")
endif()

add_subdirectory(src)

# Main executable
add_executable(TrajectoryClustering_run main.cpp)

set(GNU_OPTIONS
    "-pedantic" "-Wall" "-Wextra" "-Wcast-align" "-Wcast-qual" "-Wlogical-op"
    "-Wctor-dtor-privacy" "-Wdisabled-optimization" "-Wformat=2" "-Winit-self"
    "-Wmissing-declarations" "-Wshadow" "-Wmissing-include-dirs" "-Wnoexcept"
    "-Wold-style-cast" "-Woverloaded-virtual" "-Wredundant-decls" "-Wsign-promo"
    "-Wsign-conversion" "-Wstrict-null-sentinel" "-Wstrict-overflow=2" "-Wundef"
    "-Wswitch-default")

set(CLANG_OPTIONS
    "-Weverything" "-Wno-c++98-compat" "-Wno-covered-switch-default"
    "-Wno-float-equal" "-Wno-double-promotion" "-Wno-unused-template"
    "-Wno-c++98-compat-pedantic" "-Wno-padded" "-Wno-weak-template-vtables" 
    "-Wl,--no-as-needed" "-lprofiler -Wl" "--as-needed"
     "-fprofile-instr-generate" "fprofile-arcs" "-fcoverage-mapping")

# set(MSVC_OPTIONS "/W4")

if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    target_compile_options(TrajectoryClustering_lib PRIVATE ${GNU_OPTIONS})
    target_compile_options(TrajectoryClustering_run PRIVATE ${GNU_OPTIONS})
elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
    target_compile_options(TrajectoryClustering_lib PRIVATE ${CLANG_OPTIONS})
    target_compile_options(TrajectoryClustering_run PRIVATE ${CLANG_OPTIONS})
    target_link_options(TrajectoryClustering_lib PUBLIC "-stdlib=libstdc++")
    # Disable this if you want to use libc++
else()
    target_compile_options(TrajectoryClustering_lib PRIVATE ${MSVC_OPTIONS})
    target_compile_options(TrajectoryClustering_run PRIVATE ${MSVC_OPTIONS})
endif()

find_package(OpenMP 3.0)
if(OpenMP_CXX_FOUND)
    target_link_libraries(TrajectoryClustering_lib PUBLIC OpenMP::OpenMP_CXX)
endif()

target_link_libraries(TrajectoryClustering_run PRIVATE TrajectoryClustering_lib)
