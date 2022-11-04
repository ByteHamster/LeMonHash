cmake_minimum_required(VERSION 3.16)
project(DirectRankStoring)

if(CMAKE_BUILD_TYPE STREQUAL "Release" AND CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
    add_compile_options(-march=native)
endif()

set(BUILD_EXAMPLES OFF CACHE BOOL "Use PGM as library only")
set(BUILD_PGM_TUNER OFF CACHE BOOL "Use PGM as library only")
set(BUILD_PGM_BENCHMARK OFF CACHE BOOL "Use PGM as library only")
add_subdirectory(extlib/pgm-index EXCLUDE_FROM_ALL)

add_subdirectory(extlib/util EXCLUDE_FROM_ALL)
add_subdirectory(extlib/simpleRibbon EXCLUDE_FROM_ALL)

add_library(DirectRankStoring INTERFACE)
target_compile_features(DirectRankStoring INTERFACE cxx_std_20)
target_include_directories(DirectRankStoring INTERFACE include)
target_link_libraries(DirectRankStoring INTERFACE SimpleRibbon ByteHamsterUtil pgmindexlib)

# Do not create benchmark targets when this is a subdirectory
if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
    add_executable(IntegerDistributions benchmark/integerDistributions.cpp)
    target_compile_features(IntegerDistributions PRIVATE cxx_std_20)
    target_link_libraries(IntegerDistributions PRIVATE DirectRankStoring tlx)
endif()
