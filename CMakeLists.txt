cmake_minimum_required(VERSION 2.6.0)
project (prime-test)

set(CMAKE_BUILD_TYPE Release)

# Bring the headers
include_directories(include)

# Adding all sources
file(GLOB SOURCES "src/*.c")

add_executable (chebyshev ${SOURCES})

install(TARGETS chebyshev DESTINATION /usr/local/bin)
