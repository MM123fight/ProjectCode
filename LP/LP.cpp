//
// Created by Lu Meng on 2018/11/27.
//

#include "LPuse.h"
#include "LP.h"
int main() {
    std::string root_directory = "/Users/lumeng/Desktop/CCode/LP/data";
    const std::string all_experiments[10] = {"protein","test3","Mat0|1000_10000_0.01_0.5","Mat0|10000_60000_0.001_0.5"};
    std::string experiment = all_experiments[0] + "/SVMtoLP";
    std::string experiment = all_experiments[1];

    LP<unsigned int, double> *inst = new LP<unsigned  int, double>(root_directory, experiment);


    double rho =1;
    double tol = 1e-6;
    lp_param.sub_tol = 1e-6;
    std::string subtype = "exact";
    //subtype = "NU_ACDM";
    Print(subtype);
    Print("ADMM1");
    Print(inst->ADMM1(tol,subtype));


    delete inst;
    return 0;
}