//
// Created by Lu Meng on 2018/11/27.
//

#ifndef LP_PARAM_H
#define LP_PARAM_H

#include "../use/useheader.h"
class LP_param{
public:
    double tol = 1e-10;
    unsigned int precision = 10;
    double max_time = 100;
    //per pass
    double per_output = 1.;
    bool screenprint = true;
    double sub_tol = 1e-3;

    int inv_bound = 1000;
    //parameters for line search(SSN)
    double sigma = 0.1;//SSN tau
    double mu = 0.49;
    double delta = 0.9;
    double SSN_tau = 0.1;

    //parameters for ADMM;
    double rho = 1;
    double tau = 1;



    LP_param(){}
    virtual ~LP_param() {}
};

LP_param lp_param;
#endif //LP_PARAM_H
