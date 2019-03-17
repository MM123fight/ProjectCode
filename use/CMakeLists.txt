cmake_minimum_required(VERSION 3.10)
project(use)

set(CMAKE_CXX_STANDARD 11)



include_directories(.)

add_executable(use
        test.cpp
        Print.h
        stdafx.h Matrix.h Sample.h DataGen.h usefun.h DataRW.h)

find_package(GSL REQUIRED)
find_package(Eigen3 REQUIRED)
target_link_libraries(use GSL::gsl GSL::gslcblas Eigen3::Eigen)
target_link_libraries(use quadprog)