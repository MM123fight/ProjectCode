//
// Created by Lu Meng on 2018/12/1.
//

#ifndef LP_NU_ACDM_H
#define LP_NU_ACDM_H

#include "Param.h"

//(AT*A+scale*I)x = d
// when we update grad, we need to calculate A*delta_x, which requires AT
/*
 * x = (A*AT+scale*I)^{-1}d = scale^{-1}(I-A*(AT*A+scale*I)^{-1}AT)d
 *   = 1/scale*(d - A*(AT*A+scale*I)^{-1}ATd)
 *    (AT*A+scale*I)z = ATd
 * x = 1/scale(d - A*z)
 * (AT*AT+scale*I)x -d = 1/scale(A*AT+scale*I)(d - A*z) -d  = 1/scale(A*AT*d -A*(AT*A+scale*I)*z)
 *                = 1/scale*A*(AT*d - (AT*A+scale*I)*z) = 1/scale*A*r;
 */
template <typename L, typename D>
void NU_ACDM(std::vector<D> &x, std::vector<D> &w, D& time,const L &n, const L &m, const std::vector<L> &AT_row_ptr,
             const std::vector<L> &AT_col_idx, const std::vector<D> &AT_value, const std::vector<D> &step,
             const std::vector<D> &d, const std::vector<D> &p, const std::vector<D> &p_sum, const D &tol,
             const D &tau, const D &eta_sigma, D scale = 1, std::string Mat_type ="T", std::string cond_type = "grad_cord"){
    //time_variables
    D time_start;
    //variables
    std::vector<D> y = x;
    std::vector<D> z = x;
    std::vector<D> G(n,0);
    //G = (ATA+scale*I)x-d;
    ATx(w, n, m, AT_row_ptr, AT_col_idx, AT_value, x);
    Ax(G, n, AT_row_ptr, AT_col_idx, AT_value, w);
    for (L row = 0; row < n; ++row)
        G[row] += (scale*x[row] - d[row]);

    std::vector<L> cond_cord(n,1);
    L sum_cond_cord = 0;
    if(cond_type == "grad_cord"){
        for (L row = 0; row < n; ++row) {
            if(G[row] * G[row] <= tol){
                cond_cord[row] = 0;
            }
        }
        for (L row = 0; row < n; ++row)
            sum_cond_cord += cond_cord[row];
    }

    D grad;
    gsl_rng_env_setup();
    gsl_rng *gsl_rng_r = gsl_rng_alloc(gsl_rng_default);
    L i;
    D tmp;
    D delta_x;
    D cond_tmp;
    time_start = clock();
    while( sum_cond_cord > 0 ){
        //sample i with p
        i = RandSample(n,p_sum,gsl_rng_r) - 1;

        grad = 0;
        for (L col = AT_row_ptr[i]; col < AT_row_ptr[i+1]; ++col)
            grad += AT_value[col]*w[AT_col_idx[col]];
        grad += scale*x[i] - d[i];
        y[i] = x[i] - grad/step[i];
        z[i] = (z[i] + eta_sigma*(x[i] - 1./(scale*p[i])*grad))/(1+eta_sigma);
        tmp = tau*z[i] +(1-tau)*y[i];
        delta_x = tmp - x[i];
        x[i] = tmp;
        for (L col = AT_row_ptr[i]; col < AT_row_ptr[i+1]; ++col) {
            w[AT_col_idx[col]] += delta_x * AT_value[col];
        }

        if(cond_type == "grad_cord") {
            sum_cond_cord -= cond_cord[i];
            if (grad*grad <= tol)
                cond_cord[i] = 0;
            sum_cond_cord += cond_cord[i];
        }
    }
    time += (D)(clock() - time_start)/CLOCKS_PER_SEC;
    gsl_rng_free(gsl_rng_r);
}

template <typename L, typename D>
void NU_ACDM_param(std::vector<D> &p, std::vector<D> &p_sum, D &tau, D &eta_sigma, const L &n,
                   const std::vector<D> &step, D scale = 1, D beta = 0) {
    D alpha = (1 - beta) / 2;
    D Salpha;
    p.clear();
    p.resize(n, 0);
    p_sum.clear();
    p_sum.resize(n, 0);
    for (L row = 0; row < n; ++row)
        p[row] = pow(step[row], alpha);
    Salpha = 0;
    for (L row = 0; row < n; ++row)
        Salpha += p[row];
    for (L row = 0; row < n; ++row)
        p[row] /= Salpha;

    for (L row = 0; row < n; ++row)
        p_sum[row + 1] = p_sum[row] + p[row];

    tau = 2. / (1 + sqrt(4 * Salpha * Salpha/ scale + 1));
    eta_sigma = scale / (tau * Salpha*Salpha);
}







#endif //LP_NU_ACDM_H
