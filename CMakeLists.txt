project(vptree-draw)
cmake_minimum_required(VERSION 3.11)

add_subdirectory(src)

add_executable(poster src/poster.cpp)
add_executable(compare src/compare.cpp)

target_compile_features(poster PRIVATE cxx_std_17)
target_include_directories(poster PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/CavalierContours/include/>)

target_compile_features(compare PRIVATE cxx_std_17)
target_include_directories(compare PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/CavalierContours/include/>)
