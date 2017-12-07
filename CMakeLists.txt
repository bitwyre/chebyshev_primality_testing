cmake_minimum_required(VERSION 2.6.0)
project (chebyshev)

set(CMAKE_BUILD_TYPE Release)

SET(CMAKE_C_COMPILER /usr/local/bin/gcc-7)
SET(CMAKE_CXX_COMPILER /usr/local/bin/g++-7)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Bring the headers
include_directories(include)

# Adding all sources
file(GLOB SOURCES_C "src/*.c")
file(GLOB SOURCES_CPP "src/*.cpp")

add_executable (chebyshev ${SOURCES_C} ${SOURCES_CPP})

install(TARGETS chebyshev DESTINATION /usr/local/bin)
