//
// Created by Lu Meng on 2018/12/2.
//

#include "NU_ACDM.h"
#include "LP.h"
#include "LPuse.h"
#include "../use/Matrix.h"
#include "../use/DataRW.h"
int main() {
    std::string root_directory = "/Users/lumeng/Desktop/CCode/LP/data";
    const std::string all_experiments[10] = {"test1","Mat10|100_1000_0.1_0.5"};
    std::string experiment = all_experiments[1];

    LP<int, double> *inst = new LP<int, double>(root_directory, experiment);
    int m = inst->mi;
    int n = inst->n;
    std::cout.precision(3);
    std::vector<double> x(n,0);
    double time;
    std::vector<int> AT_row_ptr;
    std::vector<int> AT_col_idx;
    std::vector<double> AT_value;
    std::vector<double> step(n,0);
    std::vector<double> d(n,1);

    Mtranspose(m,n,inst->Ai_row_ptr,inst->Ai_col_idx,inst->Ai_value,AT_row_ptr,AT_col_idx,AT_value);
    Mdiag(n,AT_row_ptr,AT_value,step);
    for (int i = 0; i < n; ++i)
        step[i] +=1;

    std::vector<double> p(n,0);
    std::vector<double> p_sum(n+1,0);
    double tol = 1e-3;
    double tau;
    double eta_sigma;

    std::vector<double> Atimesx(m,0);

    NU_ACDM_param(p,p_sum,tau,eta_sigma,n,step);
    time  = 0;
    NU_ACDM(x,Atimesx,time, n, m, AT_row_ptr, AT_col_idx, AT_value, step, d, p, p_sum,tol, tau, eta_sigma);
    Print(time);
    Print(x);
    std::vector<double> w;
    std::vector<double> diff;
    double diff_val = 0;
    Ax(w,m,inst->Ai_row_ptr,inst->Ai_col_idx,inst->Ai_value,x);
    Axminusb(diff,n,AT_row_ptr,AT_col_idx,AT_value,w,d);
    for (int row = 0; row < n; ++row) {
        diff[row] += x[row];
        diff_val += diff[row] * diff[row];
    }
    Print("diff_val", diff_val);



    double rho =1;

    //Print(inst->ADMM1(rho,tol));

    delete inst;
    return 0;
}

