//
// Created by Lu Meng on 2018/12/2.
//

#include "NU_ACDM.h"
#include "LP.h"
#include "LPuse.h"
#include "../use/Matrix.h"
#include "../use/DataRW.h"
#include "CG.h"
#include "LS_exact.h"
#include "SSN.h"

int main() {
    std::string root_directory = "/Users/lumeng/Desktop/CCode/LP/data";
    const std::string all_experiments[10] = {"test1","Mat0|50_10000_0.1_0.5"};
    std::string experiment = all_experiments[1];

    LP<int, double> *inst = new LP<int, double>(root_directory, experiment);
    int n = inst->n;
    std::vector<double> x(n,0);
    std::vector<double> d(n);
    for (int row = 0; row < n; ++row) {
        d[row] = 1;
    }
    std::vector<double> w(inst->mi);
    for (int row = 0; row < inst->mi; ++row) {
        w[row] = -inst->bi[row];
    }
    std::vector<double> v = inst->c;
    double time = 0;

    //for (int i = 0; i < 10; ++i) {
        x.resize(n, 0);
        CG_inv(x, time, inst->mi, inst->n, inst->Ai_row_ptr, inst->Ai_col_idx, inst->Ai_value, d, 1e-15);
    //}
    Print("time CG_inv", time);
    Print(x);

    //for (int i = 0; i < 10; ++i) {
        x.resize(n,0);
        CG(x, time, inst->mi, inst->n, inst->Ai_row_ptr, inst->Ai_col_idx, inst->Ai_value, d, 1e-15);
    //}
    Print("time CG", time);
    Print(x);

    time = 0;
    for (int i = 0; i < 10; ++i) {
        x = d;
        CG(x, time, inst->mi, inst->n, inst->Ai_row_ptr, inst->Ai_col_idx, inst->Ai_value, d, 1e-15);
    }
    Print("time CG:d", time);
    Print(x);

    time = 0;
    //for (int i = 0; i < 10; ++i) {
        LS_exact(x, time, inst->mi, inst->n, inst->Ai_row_ptr, inst->Ai_col_idx, inst->Ai_value, d);
    //}
    Print("exact:", time);
    Print(x);



    //Print(inst->ADMM1(rho,tol));

    delete inst;
    return 0;
}

