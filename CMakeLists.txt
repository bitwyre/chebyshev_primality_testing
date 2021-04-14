cmake_minimum_required(VERSION 2.6.0)
project (chebyshev)

set(CMAKE_BUILD_TYPE Release)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

## Compiler flags
if(CMAKE_COMPILER_IS_GNUCXX)
    set(CMAKE_CXX_FLAGS "-O2")        ## Optimize
    set(CMAKE_CXX_FLAGS "-O3")        ## Optimize
    set(CMAKE_EXE_LINKER_FLAGS "-s")  ## Strip binary
endif()

# Bring the headers
include_directories(include)

# Adding all sources
set(SOURCES ${CMAKE_SOURCE_DIR}/src/chebyshev-polynomial.cpp
            ${CMAKE_SOURCE_DIR}/include/chebyshev-polynomial.h
            ${CMAKE_SOURCE_DIR}/include/benchmark.h
            ${CMAKE_SOURCE_DIR}/src/chebyshev-primality-test.cpp
            ${CMAKE_SOURCE_DIR}/src/primality-test-baseline.c
            ${CMAKE_SOURCE_DIR}/include/primality-test-baseline.h)

add_executable(chebyshev ${SOURCES})
target_link_libraries(chebyshev ${CMAKE_SOURCE_DIR}/libbenchmark.a pthread)