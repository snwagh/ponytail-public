cmake_minimum_required(VERSION 2.6)

project(Dolphin)

add_compile_options(-O3)
add_compile_options(-Wall)
add_compile_options(-std=c++17)
add_compile_options(-march=native)
add_compile_options(-g)

include_directories(${CMAKE_SOURCE_DIR}/external/ntl/include)
include_directories(${CMAKE_SOURCE_DIR}/external/gmp/include)
add_subdirectory(${CMAKE_SOURCE_DIR}/external/seal/native/src)

add_executable("main" "src/seal_augment.cpp" "src/main.cpp" "src/matmul.cpp")
target_link_libraries("main" seal)
target_link_libraries("main" ${CMAKE_SOURCE_DIR}/external/ntl/src/ntl.a)
target_link_libraries("main" ${CMAKE_SOURCE_DIR}/external/gmp/.libs/libgmp.a)

add_executable("test_ntt" "src/test_ntt.cpp" "src/seal_augment.cpp")
target_link_libraries("test_ntt" seal)
target_link_libraries("test_ntt" ${CMAKE_SOURCE_DIR}/external/ntl/src/ntl.a)
target_link_libraries("test_ntt" ${CMAKE_SOURCE_DIR}/external/gmp/.libs/libgmp.a)
