project(vptree-draw)
cmake_minimum_required(VERSION 3.11)

add_executable(vptree main)

target_compile_features(vptree PRIVATE cxx_std_17)
target_include_directories(vptree PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/CavalierContours/include/>)
