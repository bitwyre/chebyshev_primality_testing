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
file(GLOB SOURCES_C "src/*.c")
file(GLOB SOURCES_CPP "src/*.cpp")

add_executable (chebyshev ${SOURCES_C} ${SOURCES_CPP})

install(TARGETS chebyshev DESTINATION /usr/local/bin)
