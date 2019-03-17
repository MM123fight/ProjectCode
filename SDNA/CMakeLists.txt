cmake_minimum_required(VERSION 3.10)
project(SDNA)

set(CMAKE_CXX_STANDARD 11)

set(USE_FILES
        ../use/useheader.h)
set(LBFGSB_FILES

        ../lbfgsb/lbfgsb.c
        ../lbfgsb/linesearch.c
        ../lbfgsb/linpack.c
        ../lbfgsb/miniCBLAS.c
        ../lbfgsb/print.c
        ../lbfgsb/subalgorithms.c
        ../lbfgsb/timer.c)

set(SOURCE_FILES GenInfo.cpp ${LBFGSB_FILES} ${USE_FILES} GenInfo.h LS.h ProbParam.h Hinge.h L2LS.h Print.h)

add_executable(SDNA ${SOURCE_FILES})
find_package(GSL REQUIRED)
find_package(Eigen3 REQUIRED)
target_link_libraries(SDNA GSL::gsl GSL::gslcblas Eigen3::Eigen)
