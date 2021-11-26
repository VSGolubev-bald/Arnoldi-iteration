cmake_minimum_required(VERSION 3.20)
project(Arnoldi)

set(CMAKE_CXX_STANDARD 17)

add_executable(Arnoldi Matrix.h Matrix.cpp main.cpp)
