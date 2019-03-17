//
// Created by Lu Meng on 2018/5/10.
//

#ifndef SDNA_PROBPARAM_H
#define SDNA_PROBPARAM_H

#include "../use/useheader.h"
class ProbParam{
public:
    std::string data_dir;

    double tol = 1e-10;
    unsigned int precision = 10;
    double max_time = 100;
    //per pass
    double per_output = 1.;
    bool screenprint = true;
    double theta = 0.01;
    double eta = 0.1;
    double lambda_multiple = 1;
    double leverage_approx = 0.5;
    double leverage_prob = 0.9;


    ProbParam(){}
    virtual ~ProbParam() {}
};

ProbParam param;

#endif //SDNA_PROBPARAM_H
