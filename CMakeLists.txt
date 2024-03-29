cmake_minimum_required(VERSION 3.16)
project(LeMonHash)

if((CMAKE_BUILD_TYPE STREQUAL "Release" OR CMAKE_BUILD_TYPE STREQUAL "RelWithDebInfo")
       AND CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
    add_compile_options(-march=native)
endif()

set(TLX_INSTALL_INCLUDE_DIR tlx CACHE PATH "Workaround for TLX breaking the first cmake call")
add_subdirectory(extlib/tlx EXCLUDE_FROM_ALL)
add_subdirectory(extlib/util EXCLUDE_FROM_ALL)
add_subdirectory(extlib/simpleRibbon EXCLUDE_FROM_ALL)

add_library(PTHash INTERFACE)
target_compile_features(PTHash INTERFACE cxx_std_20)
target_include_directories(PTHash INTERFACE extlib/pthash/include)

option(SUCCINCT_USE_INTRINSICS "Use a set of intrinsics available on all x86-64 architectures" ON)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/extlib/ds2i/succinct/succinct_config.hpp.in
        ${CMAKE_CURRENT_SOURCE_DIR}/extlib/ds2i/succinct/succinct_config.hpp)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/extlib/ds2i/ds2i_config.hpp.in
        ${CMAKE_CURRENT_SOURCE_DIR}/extlib/ds2i/ds2i_config.hpp ESCAPE_QUOTES)
add_library(ds2i INTERFACE)
target_include_directories(ds2i SYSTEM INTERFACE extlib/ds2i)

add_library(LeMonHash INTERFACE)
target_compile_features(LeMonHash INTERFACE cxx_std_20)
target_include_directories(LeMonHash INTERFACE include)
target_link_libraries(LeMonHash INTERFACE SimpleRibbon ByteHamsterUtil ds2i PTHash)

# Do not create benchmark targets when this is a subdirectory
if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
    add_executable(IntegerDistributions benchmark/integerDistributions.cpp)
    target_compile_features(IntegerDistributions PRIVATE cxx_std_20)
    target_link_libraries(IntegerDistributions PRIVATE LeMonHash tlx)

    add_executable(StringDistributions benchmark/stringDistributions.cpp)
    target_compile_features(StringDistributions PRIVATE cxx_std_20)
    target_link_libraries(StringDistributions PRIVATE LeMonHash tlx)
endif()
