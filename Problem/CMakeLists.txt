cmake_minimum_required(VERSION 3.10)
project(ProbData)

set(CMAKE_CXX_STANDARD 11)

set(USE_FILES
        ../use/useheader.h)
set(SDNA_FILES
        ../SDNA/Print.h
        ../SDNA/ProbParam.h)
set(LBFGSB_FILES

        ../lbfgsb/lbfgsb.c
        ../lbfgsb/linesearch.c
        ../lbfgsb/linpack.c
        ../lbfgsb/miniCBLAS.c
        ../lbfgsb/print.c
        ../lbfgsb/subalgorithms.c
        ../lbfgsb/timer.c)

set(SOURCE_FILES LP.cpp ${LBFGSB_FILES} ${USE_FILES} ${SDNA_FILES} ProbData.h LPuse.h)

add_executable(ProbData ${SOURCE_FILES})
find_package(GSL REQUIRED)
find_package(Eigen3 REQUIRED)
target_link_libraries(ProbData GSL::gsl GSL::gslcblas Eigen3::Eigen)

