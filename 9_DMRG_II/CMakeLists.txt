cmake_minimum_required(VERSION 3.9)

#name of the project
project(TEBD VERSION 0.0.1 DESCRIPTION "this is TEBD method")

#adding all dependencies
set(SOURCE_FILES
  src/tebd.hpp
  src/main.cpp
  include/simple_linear_regression.h
)

add_executable(TEBD ${SOURCE_FILES})



#including eigen:
target_include_directories(TEBD
  PUBLIC
    src
    include
    ${CMAKE_CURRENT_SOURCE_DIR} ../8_DMRG_I
    /media/ziga/2_disk1/Downloads/eigen_build
)

find_package (Eigen3 3.3 REQUIRED NO_MODULE)
target_link_libraries(TEBD Eigen3::Eigen)
