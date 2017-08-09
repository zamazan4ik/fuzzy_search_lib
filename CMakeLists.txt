cmake_minimum_required(VERSION 3.7)
project(FuzzySearch)

set(CMAKE_CXX_STANDARD 14)

set(CMAKE_BUILD_TYPE Release)

include_directories(include/)


file(GLOB_RECURSE add_SRC
        "include/seqan/*.h"
        "include/seqan/*.cpp"
        )

set(SOURCE_FILES main.cpp Bitap.hpp Randl_fuzzy_search.hpp libflasm.h libflasm.cpp ${add_SRC})
add_executable(Fuzzyearch ${SOURCE_FILES})