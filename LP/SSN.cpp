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
    const std::string all_experiments[10] = {"test1","Mat0|100_1000_0.1_0.5"};
    std::string experiment = all_experiments[1];

    LP<int, double> *inst = new LP<int, double>(root_directory, experiment);
    int n = inst->n;
    int m = inst->mi;
    std::vector<double> x(n,0);
    std::vector<double> d(n);
    for (int row = 0; row < n; ++row) {
        d[row] = 1;
    }
    std::vector<double> w(m,0);
    for (int row = 0; row < m; ++row) {
        w[row] = -inst->bi[row];
    }

    std::vector<double> v(n,0);
    std::vector<double> AAT(m,0);
    Mdiag(m,inst->Ai_row_ptr,inst->Ai_value,AAT);
    double time = 0;
    std::string subtype= "ACDM";
    double scale_w = 1;
    double scale_v = 1;
    double tol = 1e-10;
    SSN(x,w,v,time,inst->mi, inst->n,inst->Ai_row_ptr, inst->Ai_col_idx, inst->Ai_value,AAT,inst->c,scale_w,scale_v,tol,subtype);

    delete inst;
    return 0;
}

