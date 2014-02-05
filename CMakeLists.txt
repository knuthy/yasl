cmake_minimum_required(VERSION 2.8)

project(YASL)

if($CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_BINARY_DIR})
    message(FATAL_ERROR "In-source build not allowed, please create a build directory and launch cmake in it.")
endif()

if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release")
endif()

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")

################################
# Testing
################################
# Locate GTest
set(GTEST_LIBRARY /usr/local/lib/libgtest.a)
set(GTEST_MAIN_LIBRARY /usr/local/lib/libgtest_main.a)
set(GTEST_INCLUDE_DIR /usr/local/include/)

set(GMOCK_LIBRARY /usr/local/lib/libgmock.a)
set(GMOCK_MAIN_LIBRARY /usr/local/lib/libgmock_main.a)
set(GMOCK_INCLUDE_DIR /usr/local/include/)

find_package(GTest REQUIRED)
find_package(GMock REQUIRED)
include_directories(${GTEST_INCLUDE_DIRS} ${GMOCK_INCLUDE_DIRS})

add_executable(
    runTests 
    test/test_sparse_matrix.cpp
    )
target_link_libraries(runTests
    ${GTEST_BOTH_LIBRARIES}
    ${GMOCK_BOTH_LIBRARIES}
)

enable_testing()
