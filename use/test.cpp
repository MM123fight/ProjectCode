//
// Created by Lu Meng on 2018/10/23.
//
#include <gsl/gsl_matrix.h>
#include "useheader.h"
#include "../QuadProg/src/QuadProg++.hh"
int main() {
    // The data and code directory;
    /*
    std::string root_directory = "/Users/lumeng/Desktop/Ccode/use";
    int n = 5;
    int d = 4;
    std::vector<int> A_col_idx  = {1,2,3,2,4,1,2,4,1,3,2,4};
    for (int i = 0; i < A_col_idx.size(); ++i)
        A_col_idx[i] -= 1;
    std::vector<double> A_value = {-1,2,4,1,1,1,-2,2,2,2,1,2};
    std::vector<int> A_row_ptr = {0,3,5,8,10,12};
    std::vector<int> J = {1, 2, 3};
    std::vector<int> AJ_row_ptr;
    std::vector<int> AJ_col_idx;
    std::vector<double> AJ_value;
    int sizeJ = J.size();
    subMat(AJ_row_ptr,AJ_col_idx,AJ_value,A_row_ptr,A_col_idx,A_value,sizeJ,J);
    std::vector<double> Diag;
    Diag_ATA(Diag,sizeJ,d,AJ_row_ptr,AJ_col_idx,AJ_value);


    Print("The Print.h construct the function Print() to print data or messages.\nYou can test:");
    std::string mode = "Full";
    Print("A");
    Print(n,d,A_row_ptr,A_col_idx,A_value,mode);
    Print("AJ");
    Print(sizeJ,d,AJ_row_ptr,AJ_col_idx,AJ_value);
    Print("Diag", Diag);

    Print("The Matrix.h construct the operations for sparse matrix and vector");
    std::vector<int> AT_col_idx = std::vector<int>();
    std::vector<int> AT_row_ptr = std::vector<int>();
    std::vector<double> AT_value = std::vector<double>();

    Print("Mtranspose");
    Mtranspose(n,d,A_row_ptr,A_col_idx,A_value,AT_row_ptr,AT_col_idx,AT_value);
    Print(d,n,AT_row_ptr,AT_col_idx,AT_value,mode);

    Print("Mdiag");
    std::vector<double> AL = std::vector<double>();
    Mdiag(n,A_row_ptr,A_value,AL);
    Print(AL);

    Print("Msub");
    int tau = 5;
    std::vector<int> S = {0,1,2,3,4};
    std::vector<double> Ms = std::vector<double>();
    double alpha = 0.2;
    double beta = 0.2;
    Msub(Ms,A_row_ptr,A_col_idx,A_value,AL,tau,S,beta);
    Print(tau,Ms);

    Print("WplusATdelta");
    std::vector<double> w = std::vector<double>();
    std::vector<double> delta = std::vector<double>();
    w.resize(d,0);
    delta.resize(n,1);


    Print("Ax");
    std::vector<double> x = std::vector<double>();
    x.resize(d,1);
    Print(Ax(n,A_row_ptr,A_col_idx,A_value,x));

    Print("VplusV");
    std::vector<double> y = std::vector<double>();
    y.resize(d,2);
    Print(VplusV(d,x,y));

    gsl_rng_env_setup();
    gsl_rng *gsl_rng_r = gsl_rng_alloc(gsl_rng_default);
    Print(SampleGen(10, 4, gsl_rng_r));
    gsl_rng_free(gsl_rng_r);

    Print("The DataGen generate the data, however, you should use matlab to generate");
    DataGen(10,20,0.2,root_directory);
    DataGenVector<int,double>(10,root_directory);

    Print("The DataRW read the data from files and write the data into files");
    std::vector<double> b = std::vector<double>();
    double nnzA_percent;
    MatrixRead("/Users/lumeng/Desktop/CCode/data/LS/LS5_4/LS5_4.txt",n,d,A_row_ptr,A_col_idx,A_value,b);
    Print(n,d,A_row_ptr,A_col_idx,A_value);
    Print(b);


    // min {1/2<x,Gx> + <g0,x>}
    // s.t. CE*x + ce0 = 0;
    //      CI*x + ci0 > = 0;
    Print("Use QuadProg method to solve the QP");
    std::vector<double> tmp;
    quadprogpp::Matrix<double> G, CE, CI;
    quadprogpp::Vector<double> g0, ce0, ci0, xx;
    int m, p;
    double sum = 0.0;
    char ch;
    n = 2;

    tmp = {4,-2,-2,4};
    G = VtoM_qp(n,n,tmp);
    Print(n,n,G);
    Print();

    tmp = {6,0};
    g0 = VtoV_qp(n,tmp);
    Print(n,g0);
    Print();

    m = 1;

    tmp = {1,1};
    CE = VtoM_qp(n,m,tmp);
    Print(n,m,CE);
    Print();

    tmp = {-3.0};
    ce0 = VtoV_qp(m,tmp);
    Print(m,ce0);
    Print();


    tmp = {1,0,1,0,1,1};
    p = 3;
    CI = VtoM_qp(n,p,tmp);
    Print(n,p,CI);
    Print();

    tmp = {0,0,-2};
    ci0 = VtoV_qp(p,tmp);
    Print(p,ci0);
    Print();

    xx.resize(n);

    std::cout << "f: " << solve_quadprog(G, g0, CE, ce0, CI, ci0, xx) << std::endl;
    std::cout << "x: " << xx << std::endl;

    quadprogpp::Vector<int*> k(15);
    std::vector<int> z(15);
    //int *z = new int[15];
    for (int j = 0; j<15; j++)
    {
        z[j] = j;
        k[j] = &z[j];
    }
    int size = 15;

    for (int i = 0; i < 15;i++)
    {
        std::cout << *k[i]<< " ";//因为向量容器里面都是int型的指针变量，
    }
    Print();
    z[10] = 100;
    for (int i = 0; i < 15;i++)
    {
        std::cout << *k[i]<< " ";//因为向量容器里面都是int型的指针变量，
    }
    Print();

    *k[1] = 10000;
    Print(z);
    //delete[]z;


    Print("Check LeverageA");
    n = 5;
    d = 4;
    int Bn;
    std::vector<int> B_row_ptr;
    std::vector<int> B_col_idx;
    std::vector<double> B_value;
    std::vector<double> Leverage = {0.9990,0.3437,0.9757,0.9961,0.6854};

    //1-d^{-c/3} = 1-10e-z
    // c = 3*z log(10)/log(d) we set z = 2
    //(1-epsilon)/(1+epsilon) = 10e-z = beta
    //epsilon = (1-beta)/(1+beta)
    double c = 6*log(10)/log(d);
    beta = 0.8;
    double epsilon = (1-beta)/(1+beta);
    Print("epsilon:");
    Print(epsilon);
    Print(pow(epsilon,-2));
    Print(6*log(10)*pow(epsilon,-2));
    Print(n, d,A_row_ptr,A_col_idx,A_value);
    Print("Leverage");
    LeverageA(Bn,B_row_ptr,B_col_idx,B_value,n,d,A_row_ptr,A_col_idx,A_value,Leverage,epsilon,c);
    Print(Bn, d, B_row_ptr,B_col_idx,B_value);

    x = {1,2,3,4};
    findMax<unsigned int, double>(x);
    Print(x);
    Print(findMax<unsigned int, double>(x));
    */

    Print("Check gsl_inverse:");
    std::vector<double> H = {1,0,0,2};
    std::vector<double> Hinv(4);
    gsl_matrix_view A = gsl_matrix_view_array(&H[0], 2, 2);
    gsl_matrix_view inv = gsl_matrix_view_array(&Hinv[0], 2, 2);
    gsl_matrix_inv(&A.matrix,&inv.matrix);
    Print(2,H);
    Print(2,Hinv);

    /*
    std::vector<double> prob = {0.1,0.2,0.3,0.1,0.3};
    std::vector<double> prob_sum(prob.size()+1,0);
    prob_sum[0] = 0;
    for (int row = 0; row < prob.size(); ++row) {
        prob_sum[row+1] = prob_sum[row] + prob[row];
    }
    Print(prob_sum);

    gsl_rng_env_setup();
    gsl_rng *gsl_rng_r1 = gsl_rng_alloc(gsl_rng_default);
    std::vector<int> tmp_sample(prob.size(),0);
    n = 1000;
    for (int l = 0; l < n; ++l) {
        for (int i = 0; i < prob.size(); ++i) {
            if (RandSample(prob.size(), prob_sum, gsl_rng_r1) == i+1)
                ++tmp_sample[i];
        }
    }
    for(int i = 0; i < prob.size(); ++i){
        Print((double)tmp_sample[i]/n);
    }
    gsl_rng_free(gsl_rng_r1);

    std::vector<double> v(10,0);
    for (int i1 = 0; i1 < 10; ++i1)
        v[i1] = i1;
    J = {1,4,7};
    sizeJ = J.size();
    std::vector<double> vJ;
    subVec(vJ,v,sizeJ,J);
    Print(vJ);
     */
    return 0;
}