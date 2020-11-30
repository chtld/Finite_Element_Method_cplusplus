cmake_minimum_required(VERSION 3.3)

project(FEM)

include_directories(../src/ ./ /usr/include/eigen3/)
add_executable(test_fem test.cpp)