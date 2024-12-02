cmake_minimum_required(VERSION 3.28)
project(Project)

set(CMAKE_CXX_STANDARD 17)

add_executable(Project main.cpp
        LBM.cpp
        seconds.cpp)
