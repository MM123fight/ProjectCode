//
// Created by Lu Meng on 2018/12/1.
//

#ifndef LP_SSN_H
#define LP_SSN_H

#include "LS_exact.h"
#include "CG.h"
#include "../use/useheader.h"

/* min cTx + 0.5*scale_w||w+||^2 + 0.5*scale_v||v||^2
 * w(Ax) = w + A*delta_x
 * v(x) =  v + delta_x
 * scale = scale_v/scale_w;
 * Grad/scale_w=AT*w(Ax) + scale*v(x)+c/scale_w = 0
 * Hes/scale = AT*D(w)A+scale*I;
 * subproblem: (AJTAJ +scale*I)* dir + AJTw(Ax)_J + scale*v(x) +c/scale_w= 0;
 *            d = - (AJTw(Ax)J + scale*v(x)+c/scale_w)
 *            dir = 1/scale*(I-AJT(AJ*AJT +scale*I)^{-1}AJ)d
 */

template <typename L, typename D>
void SSN(std::vector<D> &x,std::vector<D> &w, std::vector<D> &v, D &time, const L &m, const L &n,
         const std::vector<L> &A_row_ptr, const std::vector<L> &A_col_idx, const std::vector<D> &A_value,
         const std::vector<D> &AAT, const std::vector<D> &c,const D &scale_w, const D&scale_v, const D &tol,
         const std::string &subtype){

    //time variables
    D time_start;
    D scale = scale_w/scale_v;
    //variables
    L sizeJ;
    std::vector<L> J;
    std::vector<L> AJ_row_ptr;
    std::vector<L> AJ_col_idx;
    std::vector<D> AJ_value;
    std::vector<D> wJ;
    std::vector<D> grad(n,0);
    D grad_cond;
    std::vector<D> d(n,0);
    std::vector<D> step(n,0);
    std::vector<D> p(n,0);
    std::vector<D> p_sum(n+1,0);
    std::vector<D> dir;
    std::vector<D> Adir(m,0);
    std::vector<D> HJ;
    std::vector<D> AJd;
    std::vector<D> sol;
    std::vector<D> AJTsol(n,0);
    std::vector<D> ACDM_dir;
    std::vector<D> AJTACDM_dir;
    D ACDM_tol;
    D ACDM_tau;
    D eta_sigma;

    D tmp;
    D tmp1;
    D tmp2;



    //parameters for SSN
    D tau = lp_param.SSN_tau;
    D rho = lp_param.rho;
    D mu = lp_param.mu;
    D delta = lp_param.delta;
    D sigma = lp_param.sigma;
    D pow_index = 1 + tau;
    D eta_before;
    D eta;
    D line_search;

    //parameters for CG;
    D CG_tol;

    //J: index where wi is positive
    //Print("w",w);
    //Print("v",v);
    for (L row = 0; row < m; ++row) {
        if(w[row] > 0) {
            J.push_back(row);
            wJ.push_back(w[row]);
        }
    }
    //Print("wJ",wJ);
    sizeJ = J.size();
    //AJ: J cols of A
    subMat(AJ_row_ptr,AJ_col_idx,AJ_value,A_row_ptr,A_col_idx,A_value,sizeJ,J);
    time_start = clock();
    //d = - (AJTw(Ax)J + scale*v(x)+c/scale_w), grad = -scale_w*d;
    ATx(grad,sizeJ,n,AJ_row_ptr,AJ_col_idx,AJ_value,wJ);
    for (L row = 0; row < n; ++row)
        grad[row] = scale_w*grad[row] + scale_v*v[row]+c[row];
    //Print("grad",grad);
    grad_cond = l2norm_square(n,grad);

    time += (D)(clock() - time_start)/CLOCKS_PER_SEC;
    int s;
    while( grad_cond > tol ){
        time_start = clock();
        for (L row = 0; row < n; ++row)
            d[row] = - grad[row]/scale_w;
        dir.resize(n,0);
        time += (D)(clock() - time_start)/CLOCKS_PER_SEC;
        if(sizeJ == 0){
            for (int row = 0; row < n; ++row)
                dir[row] = d[row]/scale;
        }else if(subtype == "exact"){
            LS_exact(dir,time,sizeJ,n,AJ_row_ptr,AJ_col_idx,AJ_value,d,scale);
        }else if(subtype == "CG"){
            //solve linear system inexactly with CG method
            time_start = clock();
            CG_tol = 0;
            for (L row = 0; row < n; ++row)
                CG_tol += pow(grad[row],pow_index);
            CG_tol = std::min(sigma,CG_tol);
            time += (D)(clock() - time_start)/CLOCKS_PER_SEC;
            CG(dir,time, sizeJ ,n ,AJ_row_ptr, AJ_col_idx, AJ_value, d, CG_tol, scale);
        }else if(subtype == "ACDM"){
            /*
             * x = (AJT*AJ+scale*I)^{-1}d = scale^{-1}(I-AJT*(AJ*AJT+scale*I)^{-1}AJ)d
             *   = 1/scale*(d - AJT*(AJ*AJT+scale*I)^{-1}AJd)
             *(AJ*AJT+scale*I)z = AJd
             * x = 1/scale(d - AJT*z)
             * (AJT*AJ+scale*I)x -d = 1/scale(AJT*AJ+scale*I)(d - AJT*z) -d  = 1/scale(AJT*AJ*d -AJT*(AJ*AJT+scale*I)*z)
             *                = 1/scale*AJT*(AJ*d - (AJ*AJT+scale*I)*z) = 1/scale*AJT*r;
             */
            Ax(AJd,sizeJ,AJ_row_ptr,AJ_col_idx,AJ_value,d);
            for (L row = 0; row < sizeJ; ++row)
                step[row] = AAT[J[row]] + 1;
            NU_ACDM_param(p,p_sum,ACDM_tau,eta_sigma,sizeJ,step);
            ACDM_tol = 0;
            for (L row = 0; row < sizeJ; ++row)
                ACDM_tau += pow(AJd[row],pow_index);
            ACDM_tol = std::min(ACDM_tau,sigma)/sizeJ;
            NU_ACDM(ACDM_dir,AJTACDM_dir,time, sizeJ, n, AJ_row_ptr, AJ_col_idx, AJ_value, step, AJd, p, p_sum,ACDM_tol,ACDM_tau, eta_sigma);
            for (L row = 0; row < n; ++row)
                dir[row] = (d[row] - AJTACDM_dir[row])/scale;

        }else{
            Print("You don't add " + subtype + " method");
            exit(0);
        }
        //Print("dir",dir);
        time_start = clock();
        Ax(Adir,m,A_row_ptr,A_col_idx,A_value,dir);
        //Print("Adir",Adir);
        tmp = -0.5*scale_w*l2norm_square(sizeJ,wJ);
        tmp1 = 0;
        for (L row = 0; row < n; ++row)
            tmp1 += dir[row] *(c[row]+scale_v*v[row] - mu*grad[row]);
        tmp2 = 0.5*scale_v*l2norm_square(n,dir);
        eta_before = 0;
        //line search
        for (int i = 0; i < 20; ++i) {
            eta = pow(delta,i);
            for (L row = 0; row < m; ++row)
                w[row] += (eta-eta_before)*Adir[row];
            for (L row = 0; row < m; ++row) {
                if(w[row] > 0) {
                    J.push_back(row);
                    wJ.push_back(w[row]);
                }
            }
            line_search = 0.5*scale_w*l2norm_square(sizeJ,wJ) + tmp+ eta*tmp1 + pow(eta,2)*tmp2;
            if(line_search <= 0)
                break;
            eta_before = eta;

        }

        for (L row = 0; row < n; ++row) {
            v[row] += eta*dir[row];
            x[row] += eta*dir[row];
        }
        //Print("v",v);
        //Print("x",x);
        sizeJ = J.size();
        time += (D)(clock() - time_start)/CLOCKS_PER_SEC;
        subMat(AJ_row_ptr,AJ_col_idx,AJ_value,A_row_ptr,A_col_idx,A_value,sizeJ,J);
        time_start = clock();
        ATx(grad,sizeJ,n,AJ_row_ptr,AJ_col_idx,AJ_value,wJ);
        for (L row = 0; row < n; ++row)
            grad[row] = scale_w*grad[row] + scale_v*v[row]+c[row];
        grad_cond = l2norm_square(n,grad);
        time += (D)(clock() - time_start)/CLOCKS_PER_SEC;
        Print("grad_cond",grad_cond);
    }
};
#endif //LP_SSN_H
