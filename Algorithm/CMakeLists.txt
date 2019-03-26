cmake_minimum_required(VERSION 3.10)
project(Algorithm)

set(CMAKE_CXX_STANDARD 11)

set(USE_FILES
        ../use/useheader.h)
set(PROBLEM_FILES
        ../problem/ProblemHeader.h)
set(LOSS_FILES
        ../LossFun/LossFunHeader.h)
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

set(ALGO_FILES AbsAlgorithm.h AbsAlgorithmSub.h PPA.h ParamPPA.h NU_ACDM.h APPROX.h PCDM.h)
set(SOURCE_FILES PPA.cpp ${LBFGSB_FILES} ${USE_FILES} ${PROBLEM_FILES} ${SDNA_FILES} ${ALGO_FILES})

add_executable(Algorithm ${SOURCE_FILES})
find_package(GSL REQUIRED)
find_package(Eigen3 REQUIRED)
target_link_libraries(Algorithm GSL::gsl GSL::gslcblas Eigen3::Eigen)