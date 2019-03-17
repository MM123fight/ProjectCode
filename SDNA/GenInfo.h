//
// Created by Lu Meng on 2018/5/10.
//

#ifndef SDNA_GENINFO_H
#define SDNA_GENINFO_H

#include "../use/useheader.h"
#include "ProbParam.h"
#include "Print.h"
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_cblas.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>


/*!
 * consider A: n*d
 * @param n: number of rows
 * @param d: number of cols
 *
 * sparse notation for matrix A: n*d
 * @param A_row_ptr: size n+1
 * @param A_col_idx: size n*d
 * @param A_value: size n*d
 *
 * sparse notation for matrix A^T: d*n (a1 a2 ... an)
 * @param AT_row_ptr: size d+1
 * @param AT_col_idx: size d*n
 * @param AT_values: size d*n
 *
 * @param b: size n
 */


template <typename L, typename D>
class GenInfo{
public:

    L n;
    L d;

    // A: n*d
    std::vector<L> A_row_ptr = std::vector<L>();
    std::vector<L> A_col_idx = std::vector<L>();
    std::vector<D> A_value = std::vector<D>();

    std::vector<D> ATL = std::vector<D>();


    // A^T: d*n
    std::vector<L> AT_row_ptr = std::vector<L>();
    std::vector<L> AT_col_idx = std::vector<L>();
    std::vector<D> AT_value = std::vector<D>();

    std::vector<L> A_row_nnz = std::vector<L>();

    L Bn;
    std::vector<L> B_row_ptr = std::vector<L>();
    std::vector<L> B_col_idx = std::vector<L>();
    std::vector<D> B_value = std::vector<D>();

    std::vector<L> BT_row_ptr = std::vector<L>();
    std::vector<L> BT_col_idx = std::vector<L>();
    std::vector<D> BT_value = std::vector<D>();

    std::vector<D> BTL = std::vector<D>();

    // R^n
    std::vector<D> b = std::vector<D>();
    D scalevalue = 1;
    /*if bound_idx = 0, then the problem is not constrained.
     *if bound_idx > 0, then the problem is constrained.
     *We can adjust the position of variables to make the first bound_idx variables to be constrained.
     */
    std::string bound_type;
    L bound_idx = 0;
    std::vector<D> lower_bound = std::vector<D>();
    std::vector<D> upper_bound = std::vector<D>();
    std::vector<D> v = std::vector<D>();

    D lambda_min;

    std::string root_path;
    std::ofstream logFile;

    GenInfo(const std::string &root_path_name, const std::string &data_file_name, const bool &apply_leverage = false);
    virtual ~GenInfo();

    void setATbeA();
    void setATandAinfo();
    void setB(const std::string &root_path_name, const std::string &data_file_name, const bool & apply_leverage);
    void set_scale_add_L(const D &scale, const D &add);
    virtual void setLowerUpper(){}
    void setv_prepare(D &beta_sum, const L &tau);
    virtual void setv(const L &tau){}

    virtual void w_initial(std::vector<D> &w, const std::vector<D> &x)const{}
    virtual D PrimalDualGap(D &primal,D &dual, const std::vector<D> &x, const std::vector<D> &w)const{
        D gap = 0;
        return gap;
    }
    virtual void Gsub(std::vector<D> &Gs, const L &tau, const std::vector<L> &S,
                const std::vector<D> &x, const std::vector<D> &w)const{}
    virtual void w_update(std::vector<D> &w, const L &tau, const std::vector<L> &S,const std::vector<D> &delta_x)const{}


    // Set the result file to store the experiment details.
    void setResultFile(const std::string &result_path_file);
    std::vector<D> SDNA(const L &tau,const std::string &sub_type);
    void SDNA_sub_solve(std::vector<D> &h, D &time, const L &tau, const std::vector<L> &S, const L &bound_idx,
                        const std::vector<D> &Hs, const std::vector<D> &Hinv, const std::vector<D> &Gs, const std::vector<D> &x, const std::string &sub_type);
    std::vector<D> PCDM(const L &tau, const std::string &sub_type);
    std::vector<D> APPROX(const L &tau);

    std::vector<D>  ExpectMatrixInverseEinv( const std::string &root_path_name, const std::string &data_file_name);
};

template <typename L, typename D>
GenInfo<L, D>::GenInfo(const std::string &root_path_name, const std::string &data_file_name,const bool &apply_leverage) {
    root_path = root_path_name;
    std::string full_data_path = root_path_name + "/"+ data_file_name;

    Print("Loading Data...");
    MatrixRead(full_data_path,n,d,A_row_ptr,A_col_idx,A_value,b);
    Print("Finish loading data from:", full_data_path);

}

template<typename L, typename D>
GenInfo<L, D>::~GenInfo() {
    if(logFile.is_open()) {
        logFile.close();
    }
}

template<typename L, typename D>
void GenInfo<L, D>::setATbeA() {
    transform(n,d);
    AT_row_ptr = A_row_ptr;
    AT_col_idx = A_col_idx;
    AT_value = A_value;
}

template<typename L, typename D>
void GenInfo<L, D>::setATandAinfo() {
    Print("Calculate A_row_nnz..");
    A_row_nnz.clear();
    A_row_nnz.resize(n);
    for (L row = 0; row < n ; ++row)
        A_row_nnz[row] = A_row_ptr[row+1] - A_row_ptr[row];

    Print("Calculate ATA..");
    Mdiag(d,AT_row_ptr,AT_value,ATL);
}

template<typename L, typename D>
void GenInfo<L, D>::setB(const std::string &root_path_name, const std::string &data_file_name, const bool & apply_leverage) {
    if(apply_leverage == true){
        D epsilon = (1-param.leverage_approx)/(1+param.leverage_approx);
        //Print("epsilon:");
        // Print(epsilon);
        D c = -3*log(1-param.leverage_prob)/log(d);
        //Print("c:");
        //Print(c);
        std::vector<D> Leverage;
        std::string Leverage_path = root_path_name + "/leverage_" +  data_file_name;
        VectorRead(Leverage_path,n,Leverage);
        Print("Calculate Matrix Approximation");
        LeverageA(Bn,B_row_ptr, B_col_idx,B_value,n,d,A_row_ptr,A_col_idx,A_value,Leverage,epsilon,c);
        Print("Finish calculating");
        Print("Calculate BT..");
        Mtranspose(Bn,d,B_row_ptr,B_col_idx,B_value,BT_row_ptr,BT_col_idx,BT_value);
        Print("Finish calculating.");
        Print("Calculate diagonal entries of BTB");
        Mdiag(d,BT_row_ptr,BT_value,BTL);
        Print("Finish calculating.");

    }else{
        Bn = n;
        B_row_ptr = A_row_ptr;
        B_col_idx = A_col_idx;
        B_value = A_value;
        BT_row_ptr= AT_row_ptr;
        BT_col_idx = AT_col_idx;
        BT_value = AT_value;
        BTL = ATL;
    }

}

template<typename L, typename D>
void GenInfo<L, D>::setv_prepare(D &beta_sum,const L &tau){
    D eta = (tau-1)/(D)std::max<L>(1,d-1);
    std::vector<D> beta(n);
    beta_sum = 0;
    for (L row = 0; row < n; ++row) {
        beta[row] = 1 + eta * (A_row_nnz[row] - 1);
        beta_sum += beta[row];
    }
    v.clear();
    v.resize(d,0);
    for (L row = 0; row < d; ++row) {
        for (L col = AT_row_ptr[row]; col < AT_row_ptr[row + 1]; ++col)
            v[row] += beta[AT_col_idx[col]] * AT_value[col] * AT_value[col];
    }
}


template<typename L, typename D>
void GenInfo<L, D>::setResultFile(const std::string &result_path_file) {
    std::string full_result_path = root_path + "/result/" + result_path_file;
    logFile.open(full_result_path.c_str());
    if (logFile.fail()) {
        Print("!!! Cannot open experiment result file: ",full_result_path);
        exit(0);
    }
    Print("Set experiment result file: ",full_result_path);
}

template<typename L, typename D>
std::vector<D> GenInfo<L,D>::SDNA(const L &tau,const std::string &sub_type){
    ReportError(GenInfo<L,D>::logFile,GenInfo<L,D>::b);
    D primal;
    D dual;
    D gap;
    std::vector<D> x(d, 0);
    std::vector<D> w;
    w_initial(w,x);
    gap = PrimalDualGap(primal, dual, x, w);
    PrintInital(tau,gap,primal,dual,param.precision,logFile);

    D elapsedTime = 0;
    D start;
    D elapsedTimeSampling = 0;
    D elapsedTimeComputingHessian = 0;
    D elapsedTimeComputingGradient = 0;
    D elapsedTimeSubSolve = 0;
    D elapsedTime_x_w_update= 0;
    D elapsedTimeSub = 0;

    std::vector<L> S(tau);
    std::vector<D> Hs(tau*tau);
    std::vector<D> Gs(tau);
    std::vector<D> h(tau);

    gsl_rng_env_setup();
    gsl_rng *gsl_rng_r = gsl_rng_alloc(gsl_rng_default);
    gsl_vector *T = gsl_vector_alloc(tau);
    gsl_permutation * p = gsl_permutation_alloc(tau);
    /*This function allocates memory for a new permutation of size n.
     * The permutation is not initialized and its elements are undefined.
     */

    std::vector<D> Hinv(tau*tau);
    gsl_matrix_view H = gsl_matrix_view_array(&Hs[0], tau, tau);
    gsl_matrix_view inv = gsl_matrix_view_array(&Hinv[0], tau, tau);
    gsl_vector_view bb = gsl_vector_view_array(&Gs[0], tau);



    L itertau = 0;
    D tmp;
    while(gap > param.tol && elapsedTime <= param.max_time){
    //for (int l = 0; l < 3; ++l) {
        itertau += tau;

        start = clock();
        S = SampleGen(d, tau, gsl_rng_r);
        //Print(S);
        elapsedTimeSampling += (D)(clock() - start)/CLOCKS_PER_SEC;

        start = clock();
        Msub(Hs, BT_row_ptr, BT_col_idx, BT_value, BTL, tau, S, scalevalue);
        //Print(Hs);
        elapsedTimeComputingHessian += (D)(clock() - start)/CLOCKS_PER_SEC;
        //gsl_matrix_fprintf (stdout, &mm.matrix, "%g");

        start = clock();
        Gsub(Gs,tau,S,x,w);
        //Print(Gs);
        elapsedTimeComputingGradient += (D)(clock() - start)/CLOCKS_PER_SEC;
        //gsl_vector_fprintf (stdout, &bb.vector, "%g");

        // compute T = Qdata^{-1}GS, update xS = xS - T
        int s;
        if (sub_type == "exact") {
            start = clock();
            gsl_linalg_LU_decomp(&H.matrix, p, &s);
            gsl_linalg_LU_solve(&H.matrix, p, &bb.vector, T);
            elapsedTimeSubSolve += (D) (clock() - start) / CLOCKS_PER_SEC;

            for (L i = 0; i < tau; ++i) {
                h[i] = -T->data[i];
            }
        } else if (sub_type == "GM") {
            start = clock();
            gsl_matrix_inv(&H.matrix, &inv.matrix);
            elapsedTimeComputingHessian += (D) (clock() - start) / CLOCKS_PER_SEC;
            SDNA_sub_solve(h, elapsedTimeSubSolve, tau, S, bound_idx, Hs,Hinv,Gs, x, sub_type);
            //Print(h);
        }

        start = clock();
        for (L i = 0; i < tau; ++i)
            x[S[i]] += h[i];
        //Print(x);
        w_update(w,tau,S,h);
        //Print(w);
        elapsedTime_x_w_update += (D)(clock() - start)/CLOCKS_PER_SEC;

        elapsedTimeSub =  elapsedTimeComputingGradient + elapsedTimeSubSolve + elapsedTime_x_w_update;
        elapsedTime = elapsedTimeSampling + elapsedTimeComputingHessian+ elapsedTimeSub;
        //=========================================================================================
        // Output the data
        if (itertau % (int)(param.per_output*GenInfo<L,D>::d) == 0) {
            gap = PrimalDualGap(primal, dual,x,w);
            PrintData(itertau/GenInfo<L,D>::d, param.precision, gap, primal, dual,
                      GenInfo<L,D>::logFile, elapsedTime, elapsedTimeSub, param.screenprint);
        }

    }
    PrintFinal(elapsedTime, elapsedTimeSampling, elapsedTimeComputingHessian, elapsedTimeSub, itertau);
    GenInfo<L,D>::logFile.close();
    gsl_rng_free(gsl_rng_r);
    gsl_permutation_free(p);
    gsl_vector_free(T);
    return x;
}

template<typename L, typename D>
void GenInfo<L, D>::SDNA_sub_solve(std::vector<D> &h, D &time, const L &tau, const std::vector<L> &S, const L &bound_idx,
                                   const std::vector<D> &Hs, const std::vector<D> &Hinv, const std::vector<D> &Gs,
                                   const std::vector<D> &x, const std::string &sub_type) {
    std::vector<D> g = Gs;
    D left = 1;
    D right = 0;
    D theta = param.theta;
    D start;
    L tmp;
    L iter;
    std::vector<L> Ss = std::vector<L>();
    std::vector<D> r(tau);
    std::vector<D> Hinvr(tau);
    std::vector<D> Hh(tau);
    for (L i = 0; i < tau; ++i) {
        if(S[i] < bound_idx)
            Ss.push_back(i);
    }
    L size = Ss.size();
    L per_iter = (int) log(8*tau*tau/param.theta);
    start = clock();
    while (left > theta* right) {
        for (L row = 0; row < tau; ++row)
            h[row] -= g[row] / (tau*Hs[row]);
        for (L row = 0; row < size; ++row) {
            tmp = Ss[row];
            if (h[tmp] + x[S[tmp]] < lower_bound[S[tmp]])
                h[tmp] = lower_bound[S[tmp]] - x[S[tmp]];
            else if (h[tmp] + x[S[tmp]] > upper_bound[S[tmp]])
                h[tmp] = upper_bound[S[tmp]] - x[S[tmp]];
        }
        for (L row = 0; row < tau; ++row) {
            Hh[row] = 0;
            for (L col = 0; col < tau; ++col) {
                Hh[row] += Hs[row*tau+col] * h[col];
            }
            g[row] = Hh[row] + Gs[row];
        }

        ++iter;
        if(iter%per_iter == 0) {
            right = 0;
            for (L row = 0; row < tau; ++row) {
                right -= h[row] * (0.5 * Hh[row] + Gs[row]);
            }

            for (L row = 0; row < tau; ++row)
                r[row] = g[row];
            for (L row = 0; row < size; ++row) {
                tmp = Ss[row];
                if (g[tmp] < -1)
                    r[tmp] = g[tmp] + 1;
                else if (g[tmp] < 0) {
                    r[tmp] = 0;
                }
            }

            for (L row = 0; row < tau; ++row) {
                Hinvr[row] = 0;
                for (L col = 0; col < tau; ++col) {
                    Hinvr[row] += Hinv[row * tau + col] * r[col];
                }
            }
            left = 0.5 * VtimesV(tau, r, Hinvr);
        }
    }
    time += (D)(clock() - start)/CLOCKS_PER_SEC;

}

/*
template<typename L, typename D>
void GenInfo<L, D>::SDNA_sub_solve(std::vector<D> &h, D &time, const L &tau, const std::vector<L> &S, const L &bound_idx,
                                   const std::vector<D> &Hs,const std::vector<D> &Hinv, const std::vector<D> &Gs, const std::vector<D> &x, const std::string &sub_type){
    std::vector<D> g = Gs;
    D left = 1;
    D right = 0;
    D step = 0;
    L tmp;
    D theta = param.theta;
    L size;
    D start;
    //std::vector<D> HinvG(tau);
    std::vector<D> subgrad(tau,0);
    std::vector<D> r(tau);
    std::vector<D> h_before(tau,0);
    std::vector<D> y(tau,0);
    D rho_before = 1;
    D rho_after;
    D eta;
    D beta;
    D rho_square;
    std::vector<D> Hy(tau);
    L iter = 0;

    std::vector<D> Hinvr(tau,0);

    //HinvG = Ax(tau,Hinv,Gs);
    for (L row = 0; row < tau; ++row)
        step += Hs[row*tau+row];
    //Print("step");
   // Print(step);
    L per_iter = (int) log(8*step*step/(param.theta*lambda_min*lambda_min));
    //Print("per_iter");
    //Print(per_iter);

    std::vector<L> Ss = std::vector<L>();
    for (L i = 0; i < tau; ++i) {
        if(S[i] < bound_idx)
            Ss.push_back(i);
    }
    size = Ss.size();
    per_iter = 100;
    start = clock();
    while (left > theta* right) {
        for (L row = 0; row < tau; ++row) {
            h[row] = y[row]-g[row] / step;
        }
        for (L row = 0; row < size; ++row) {
            tmp = Ss[row];
            if (h[tmp] + x[S[tmp]] < lower_bound[S[tmp]])
                h[tmp] = lower_bound[S[tmp]] - x[S[tmp]];
            else if (h[tmp] + x[S[tmp]] > upper_bound[S[tmp]])
                h[tmp] = upper_bound[S[tmp]] - x[S[tmp]];
        }
        rho_square = rho_before*rho_before;
        eta = lambda_min/step - rho_square;
        rho_after = 0.5*(eta+sqrt(eta*eta+4*rho_square));
        beta = rho_before*(1-rho_before)/(rho_square+rho_after);

        for (L row = 0; row < tau; ++row) {
            Hy[row] = 0;
            for (L col = 0; col < tau; ++col) {
                Hy[row] += Hs[row * tau + col] * y[col];
            }
            g[row] = Hy[row] + Gs[row];
            y[row] = h[row] + beta * (h[row] - h_before[row]);
            h_before[row] = h[row];
        }
        ++iter;
        if(iter%per_iter == 0) {
            right = 0;
            for (L row = 0; row < tau; ++row) {
                right -= h[row]*(0.5*g[row] + Gs[row]);
            }

            for (L row = 0; row < tau; ++row)
                r[row] = g[row];
            for (L row = 0; row < size; ++row) {
                tmp = Ss[row];
                if (g[tmp] < -1)
                    r[tmp] = g[tmp] +1;
                else if(g[tmp] < 0){
                    r[tmp] = 0;
                }
            }
            for (L row = 0; row < tau; ++row) {
                Hinvr[row] = 0;
                for (L col = 0; col < tau; ++col) {
                   Hinvr[row] += Hinv[row*tau+col]*r[col];
                }
            }
            left = 0.5 * VtimesV(tau, r, Hinvr);
            //Print("right");
            //Print(right);
            //Print("left");
            //Print(left);
        }
    }
    time += (D)(clock() - start)/CLOCKS_PER_SEC;
}
*/

template<typename L, typename D>
std::vector<D> GenInfo<L,D>::APPROX(const L &tau) {
    setv(tau);
    D primal;
    D dual;
    D gap;

    D theta_before= tau/(D)d;
    D thetasquare_before = theta_before * theta_before ;

    std::vector<D> z(d,0);
    std::vector<D> u(d,0);
    std::vector<D> w_z;
    std::vector<D> w_u(n,0);
    w_initial(w_z,z);
    std::vector<D> x = VplusV(d,AlphatimesV(d,thetasquare_before,u),z);
    std::vector<D> w_x;
    w_initial(w_x,x);
    gap = PrimalDualGap(primal, dual, x, w_x);
    PrintInital(tau,gap,primal,dual,param.precision,logFile);

    D elapsedTime = 0;
    D start;
    D elapsedTimeSampling = 0;
    D elapsedTimeComputingGradient = 0;
    D elapsedTimeSubSolve = 0;
    D elapsedTimeSub = 0;


    std::vector<L> S(tau);
    std::vector<D> Gs(tau);
    std::vector<D> Gs_z(tau);
    std::vector<D> Gs_u(tau);
    std::vector<D> h_z(tau);
    std::vector<D> h_u(tau);


    gsl_rng_env_setup();
    gsl_rng *gsl_rng_r = gsl_rng_alloc(gsl_rng_default);

    L itertau = 0;
    D tmp;
    D eta;
    D theta_after;
   while(gap > param.tol && elapsedTime <= param.max_time){
    //for (int j = 0; j < 100; ++j) {
        itertau += tau;

        start = clock();
        S = SampleGen(d, tau, gsl_rng_r);
        elapsedTimeSampling += (D)(clock() - start)/CLOCKS_PER_SEC;

        start = clock();
        Gsub(Gs_u,tau,S,u,w_u);
        Gsub(Gs_z,tau,S,z,w_z);
        for (L row = 0; row < tau; ++row)
            Gs[row] = thetasquare_before * Gs_u[row] + Gs_z[row];
        elapsedTimeComputingGradient += (D)(clock() - start)/CLOCKS_PER_SEC;

        start = clock();
        tmp = (theta_before*d)/(D)tau;
        thetasquare_before = theta_before*theta_before;
        theta_after = 0.5*(sqrt(thetasquare_before*thetasquare_before + 4*thetasquare_before)-thetasquare_before);
        for (L i = 0; i < tau; ++i) {
            h_z[i] = -Gs[i]/(tmp*v[S[i]]);
            if(S[i] < bound_idx){
                if(h_z[i] + z[S[i]] < lower_bound[S[i]]) {
                    h_z[i] = lower_bound[S[i]] - z[S[i]];
                }
                else if(h_z[i] + z[S[i]] > upper_bound[S[i]]){
                    h_z[i] = upper_bound[S[i]] - z[S[i]];
                }
            }
            z[S[i]] += h_z[i];
            eta = (1-tmp)/thetasquare_before;
            h_u[i] = -eta*h_z[i];
            u[S[i]] += h_u[i];
        }
        w_update(w_z,tau,S,h_z);
        w_update(w_u,tau,S,h_u);
        theta_before = theta_after;
        elapsedTimeSubSolve += (D)(clock() - start)/CLOCKS_PER_SEC;
        //Only for calculating primal and dual gap, the time is not calculated;

        elapsedTimeSub =  elapsedTimeComputingGradient + elapsedTimeSubSolve;
        elapsedTime = elapsedTimeSampling + elapsedTimeSub;
        //=========================================================================================
        // Output the data
        if (itertau % (int)(param.per_output*GenInfo<L,D>::d) == 0) {
            x = VplusV(d,AlphatimesV(d,thetasquare_before,u),z);
            w_x = VplusV(d,AlphatimesV(d,thetasquare_before,w_u),w_z);
            w_initial(w_x,x);
            gap = PrimalDualGap(primal, dual,x,w_x);
            PrintData(itertau/GenInfo<L,D>::d, param.precision, gap, primal, dual,
                      GenInfo<L,D>::logFile, elapsedTime,elapsedTimeSub, param.screenprint);
        }

    }
    PrintFinal(elapsedTime, elapsedTimeSampling, elapsedTimeSub, itertau);
    GenInfo<L,D>::logFile.close();
    gsl_rng_free(gsl_rng_r);
    return x;
}


template <typename L, typename D>
std::vector<D> GenInfo<L,D>::PCDM(const L &tau, const std::string &sub_type){
    setv(tau);
    ReportError(GenInfo<L,D>::logFile,GenInfo<L,D>::b);
    L omega = findMax<L,L>(A_row_nnz);
    D beta= 1+(omega-1)*(tau-1)/(D)std::max<L>(1,d-1);
    D primal;
    D dual;
    D gap;
    std::vector<D> x(d, 0);
    std::vector<D> w;
    w_initial(w,x);
    gap = PrimalDualGap(primal, dual, x, w);
    PrintInital(tau,gap,primal,dual,param.precision,logFile);

    D elapsedTime = 0;
    D start;
    D elapsedTimeSampling = 0;
    D elapsedTimeComputingGradient = 0;
    D elapsedTimeSubSolve = 0;
    D elapsedTimeSub = 0;

    std::vector<L> S(tau);
    std::vector<D> Gs(tau);
    std::vector<D> h(tau);

    gsl_rng_env_setup();
    gsl_rng *gsl_rng_r = gsl_rng_alloc(gsl_rng_default);

    L itertau = 0;
    while(gap > param.tol && elapsedTime <= param.max_time){
        //for (int l = 0; l < 3; ++l) {
        itertau += tau;

        start = clock();
        S = SampleGen(d, tau, gsl_rng_r);
        //Print(S);
        elapsedTimeSampling += (D)(clock() - start)/CLOCKS_PER_SEC;


        start = clock();
        Gsub(Gs,tau,S,x,w);
        //Print(Gs);
        elapsedTimeComputingGradient += (D)(clock() - start)/CLOCKS_PER_SEC;

        start = clock();
        for (L i = 0; i < tau; ++i) {
            if(sub_type == "newPCDM"){
                h[i] = -Gs[i]/v[S[i]];
            }else if(sub_type == "PCDM"){
                h[i] = -Gs[i]/(beta*ATL[S[i]]);
            }else{
                Print("No such method");
            }
            //When it is constrained problem
            if(S[i] < bound_idx){
                if(h[i] + x[S[i]] < lower_bound[S[i]]) {
                    h[i] = lower_bound[S[i]] - x[S[i]];
                }
                else if(h[i] + x[S[i]] > upper_bound[S[i]]){
                    h[i] = upper_bound[S[i]] - x[S[i]];
                }
            }
            x[S[i]] += h[i];
        }

        w_update(w,tau,S,h);
        elapsedTimeSubSolve += (D)(clock() - start)/CLOCKS_PER_SEC;

        elapsedTimeSub =  elapsedTimeComputingGradient  + elapsedTimeSubSolve;
        elapsedTime = elapsedTimeSampling +  elapsedTimeSub;
        //=========================================================================================
        // Output the data


        if (itertau % (int)(param.per_output*GenInfo<L,D>::d) == 0) {
            gap = PrimalDualGap(primal, dual,x,w);
            PrintData(itertau/GenInfo<L,D>::d, param.precision, gap, primal, dual,
                      GenInfo<L,D>::logFile, elapsedTime, elapsedTimeSub, param.screenprint);
        }


    }
    PrintFinal(elapsedTime, elapsedTimeSampling, elapsedTimeSub, itertau);
    GenInfo<L,D>::logFile.close();
    gsl_rng_free(gsl_rng_r);
    return x;
}

template<typename L, typename D>
void GenInfo<L, D>::set_scale_add_L(const D &scale, const D &add) {
    for (L row = 0; row < d; ++row) {
        ATL[row] *= scale;
        ATL[row] += add;
    }
    for (L row = 0; row < d; ++row) {
        BTL[row] *= scale;
        BTL[row] += add;
    }
}

template<typename L, typename D>
std::vector<D> GenInfo<L, D>::ExpectMatrixInverseEinv( const std::string &root_path_name, const std::string &data_file_name) {
    Print("Calculate AT..");
    Mtranspose(n,d,A_row_ptr,A_col_idx,A_value,AT_row_ptr,AT_col_idx,AT_value);
    setATandAinfo();
    L times = 10000;

    D epsilon = (1 - param.leverage_approx) / (1 + param.leverage_approx);
    D c = -3 * log(1 - param.leverage_prob) / log(d);

    std::vector<D> test(n*n);
    std::vector<D> Einv(n*n);
    std::vector<D> AAT;
    std::vector<D> AL;
    Mdiag(n,A_row_ptr,A_value,AL);
    for (L row = 0; row < n; ++row) {
        AL[row] += 1;
    }
    std::vector<D> Leverage;
    std::string Leverage_path = root_path_name + "/leverage_" + data_file_name;
    VectorRead(Leverage_path, n, Leverage);
    L tau;
    std::vector<L> S;
    std::vector<D> Hs(tau*tau);
    std::vector<D> Hinv(tau*tau);
    gsl_matrix_view H = gsl_matrix_view_array(&Hs[0], tau, tau);
    gsl_matrix_view inv = gsl_matrix_view_array(&Hinv[0], tau, tau);

    for (int row = 0; row < n; ++row)
        S.push_back(row);
    Msub(AAT, A_row_ptr, A_col_idx, A_value, AL, n,S);


    for ( L i = 0; i < times; ++i) {
        leverage_sample(S, n, d, Leverage, epsilon, c);
        tau = S.size();
        Print(tau);
        for (L row = 0; row < tau; ++row) {
            for (L col = 0; col < tau; ++col) {
                Hs[row*tau + col] = AAT[S[row]*n+ S[col]];
            }
        }
        gsl_matrix_inv(&H.matrix, &inv.matrix);
        for (L row = 0; row < tau; ++row) {
            for (L col = 0; col < tau; ++col) {
                Einv[S[row]*n+ S[col]]+= Hinv[row*tau + col];
            }
        }
    }

    for (L row = 0; row < n; ++row) {
        for (L col = 0; col < n; ++col) {
            Einv[row*n+ col]/= times;
        }
    }
    for (L row = 0; row < n; ++row) {
        for (L col = 0; col < n; ++col) {
            Einv[row*n+ col]/= times;
        }
    }

    for (L row = 0; row < n; ++row) {
        for (L col = row; col < n; ++col) {
            test[row*n+col] = 0;
            for (L i = 0; i < n; ++i) {
                test[row*n+col] += AAT[row*n+i]*Einv[col*n+i];
            }
            test[col*n+row] = test[row*n+col];
        }
    }

    return test;
}


#endif //SDNA_GENINFO_H
