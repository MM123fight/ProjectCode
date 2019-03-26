//
// Created by Lu Meng on 2019/2/20.
//

#ifndef LOSSFUN_PLUS_SQUARED_LOSS_H
#define LOSSFUN_PLUS_SQUARED_LOSS_H

#include "AbsLoss.h"

/*
 * 0.5||[Ai*x -bi]_{+}||^2 + 0.5||Ae*x -be||^2  +0.5*alpha*||x-d||^2
 * = sum_{i=1}^m {0.5 ||a_ix -b_i||^2 + 0.5*alpha/m*||x-d||^2}
 *
 * wi = Ai*x - bi, we = Ae*x-be
 */

template<typename L, typename D>
class plus_squared_loss:public AbsLoss<L, D>{
private:
    std::vector<L> full_block;
public:

    plus_squared_loss(const D& alpha_value, const std::vector<D>& b_value, const std::vector<D>& d_value):
            AbsLoss<L, D>(alpha_value,b_value,d_value){
        L n = d_value.size();
        full_block = std::vector<L>(n);
        for (L row = 0; row < n; ++row)
            full_block[row] = row;
    }
    ~plus_squared_loss(){}
    //grad(x) = AiT[Ai*x-bi]_{+} +AeT(Ae*x-be) + alpha*(x-d)
    void cord_grad_update(std::vector<D>& cord_grad, const ProbData<L, D>* const data, const std::vector<D>& w,
                          const std::vector<D>& x, const L& blocksize, const std::vector<L>& cord);
    void grad_compute(std::vector<D>& grad, const ProbData<L, D>* const data,const std::vector<D>& w,
                      const std::vector<D>& x);
    void cord_grad_update(D& cord_grad, const ProbData<L, D>* const data, const D& coeff_1, const D& coeff_2,
                          const std::vector<D>& w_y_bar, const std::vector<D>& y_bar, const std::vector<D>& w_z_bar,
                          const std::vector<D>& z_bar, const L& cord);
    void grad_compute(std::vector<D> &grad, const ProbData<L, D>* const data, const D &tau,
                      const std::vector<D> &M, const std::vector<D> &w_y_bar, const std::vector<D> &y_bar,
                      const std::vector<D> &w_z_bar, const std::vector<D> &z_bar);

    /* w_y = A*y, w_z = A*z-b;
     * grad(theta_y*y + z) = A[theta_y *w_y + w_z]_{mi+} + alpha*(theta_y*y + z-d);
     */
    void cord_grad_sum_update(std::vector<D>& cord_grad_sum, const ProbData<L, D>* const data, const std::vector<D>& w_y,
                              const std::vector<D>& y, const std::vector<D>& w_z, const std::vector<D>& z,
                              const D& theta_y, const L& blocksize, const std::vector<L>& cord);
    void grad_sum_compute(std::vector<D> &grad_sum,const ProbData<L, D>* const data, const std::vector<D> &w_y, const std::vector<D> &y,
                          const std::vector<D> &w_z, const std::vector<D> &z, const D &theta_y);
    void fun_value_compute(D& fun_value, const ProbData<L, D>* const data, const std::vector<D>& w, const std::vector<D>& x);
    void v_set(std::vector<D>& v, const ProbData<L, D>* const data, const L& blocksize);
    void Lip_set(std::vector<D>& Lip, const ProbData<L, D>* const data);
};

//Update for cord for PCDM and APROX: w = Ax -b
template<typename L, typename D>
void plus_squared_loss<L, D>::cord_grad_update(std::vector<D>& cord_grad, const ProbData<L, D>* const data, const std::vector<D>& w,
                                               const std::vector<D>& x, const L& blocksize, const std::vector<L>& cord){

    L tmp_cord,tmp_idx;
    D tmp;

    for (L row = 0; row < blocksize; ++row) {
        tmp_cord = cord[row];
        cord_grad[row] = AbsLoss<L, D>::alpha * (x[tmp_cord] - AbsLoss<L, D>::d[tmp_cord]);
        for (L col = data->AT_row_ptr[tmp_cord]; col < data->AT_row_ptr[tmp_cord + 1]; ++col) {
            tmp_idx = data->AT_col_idx[col];
            tmp = w[tmp_idx];
            if ((tmp_idx < data->mplus)&&(tmp < 0))
                tmp = 0;
            cord_grad[row] += data->AT_value[col] * tmp;

        }
    }

}

template<typename L, typename D>
void plus_squared_loss<L, D>::grad_compute(std::vector<D>& grad, const ProbData<L, D>* const data,
                                           const std::vector<D>& w, const std::vector<D>& x){
    cord_grad_update(grad,data,w,x,data->n,full_block);
}

//Update for cord for NU_ACDM: w_y = Ay
//grad = Ai^T*[Ay-b]_{+} + alpha * (xi -di)
/* w_y = coeff_1 * w_y_bar + coeff_2 * w_z_bar;
 * x = coeff_1_* y_bar + coeff_2 * z_bar;
 */
template<typename L, typename D>
void plus_squared_loss<L, D>::cord_grad_update(D& cord_grad, const ProbData<L, D>* const data, const D& coeff_1, const D& coeff_2,
                                               const std::vector<D>& w_y_bar, const std::vector<D>& y_bar, const std::vector<D>& w_z_bar,
                                               const std::vector<D>& z_bar, const L& cord){
    L tmp_idx;
    D tmp;
    cord_grad = AbsLoss<L, D>::alpha * (coeff_1 * y_bar[cord] + coeff_2 * z_bar[cord]- AbsLoss<L, D>::d[cord]);
    for (L col = data->AT_row_ptr[cord]; col < data->AT_row_ptr[cord + 1]; ++col) {
        tmp_idx = data->AT_col_idx[col];
        tmp = coeff_1 * w_y_bar[tmp_idx] + coeff_1 * w_z_bar[tmp_idx] - AbsLoss<L, D>::b[tmp_idx];
        if ((tmp_idx < data->mplus) && (tmp < 0))
            tmp = 0;
        cord_grad += data->AT_value[col] * tmp;

    }

}

//grad_compute for NU_ACDM
//coeff_1 = 1/(1-tau), coeff_2 = -tau/(1-tau)
template<typename L, typename D>
void plus_squared_loss<L, D>::grad_compute(std::vector<D> &grad, const ProbData<L, D>* const data, const D &tau,
                                           const std::vector<D> &M, const std::vector<D> &w_y_bar, const std::vector<D> &y_bar,
                                           const std::vector<D> &w_z_bar, const std::vector<D> &z_bar) {
    D coeff_1 = 1./(1-tau);
    D coeff_2 = 1. - coeff_1;
    D tmp1 = coeff_1 * M[0] + coeff_2 * M[2];
    D tmp2 = coeff_1 * M[1] + coeff_2 * M[3];
    for (L row = 0; row < data->n; ++row) {
        cord_grad_update(grad[row],data,tmp1, tmp2, w_y_bar,y_bar, w_z_bar, z_bar,row);
    }

}



template<typename L, typename D>
void plus_squared_loss<L, D>::cord_grad_sum_update(std::vector<D>& cord_grad_sum, const ProbData<L, D>* const data,
                                                   const std::vector<D>& w_y, const std::vector<D>& y,
                                                   const std::vector<D>& w_z, const std::vector<D>& z,
                                                   const D& theta_y, const L& blocksize, const std::vector<L>& cord) {
    L tmp_cord,tmp_idx;
    D tmp;
    for (L row = 0; row < blocksize; ++row) {
        tmp_cord = cord[row];
        cord_grad_sum[row] = AbsLoss<L, D>::alpha * (theta_y*y[tmp_cord] + z[tmp_cord] - AbsLoss<L, D>::d[tmp_cord]);
        for (L col = data->AT_row_ptr[tmp_cord]; col < data->AT_row_ptr[tmp_cord + 1]; ++col) {
            tmp_idx = data->AT_col_idx[col];
            tmp = theta_y * w_y[tmp_idx] + w_z[tmp_idx];
            if ((tmp_idx < data->mplus)&&(tmp < 0))
                tmp = 0;
            cord_grad_sum[row] += data->AT_value[col] * tmp;

        }
    }

}

template<typename L, typename D>
void plus_squared_loss<L, D>::grad_sum_compute(std::vector<D> &grad_sum,const ProbData<L, D>* const data, const std::vector<D> &w_y, const std::vector<D> &y,
                                          const std::vector<D> &w_z, const std::vector<D> &z, const D &theta_y) {

    L n = data->n;
    L m = data->m;
    std::vector<D> sum(n);
    std::vector<D> sum_w(m);
    for (L row = 0; row < n; ++row)
        sum[row] = theta_y * y[row] + z[row];
    for (L row = 0; row < m; ++row)
        sum_w[row] = theta_y * w_y[row] + w_z[row];
    grad_compute(grad_sum,data,sum_w,sum);
}

template<typename L, typename D>
void plus_squared_loss<L, D>::fun_value_compute(D& fun_value,const ProbData<L, D> *const data, const std::vector<D>& w, const std::vector<D>& x){
    fun_value = 0;
    D tmp;
    for (L row = 0; row < data->n; ++row) {
        tmp = x[row] - AbsLoss<L, D>::d[row];
        fun_value += tmp * tmp;
    }
    fun_value *= AbsLoss<L, D>::alpha;
    for (L row = 0; row < data->m; ++row){
        tmp = w[row];
        if( (tmp > 0)||(row >= data->mplus)) {
            fun_value += tmp * tmp;
        }
    }
    fun_value *= 0.5;
};


template<typename L, typename D>
void plus_squared_loss<L, D>::v_set(std::vector<D> &v, const ProbData<L, D> *const data, const L &blocksize) {
    AbsLoss<L, D>::v_initial_set(v, data, blocksize);
    for (L row = 0; row < data->n; ++row)
        v[row] += AbsLoss<L,D>::alpha;
}

template<typename L, typename D>
void plus_squared_loss<L, D>::Lip_set(std::vector<D> &Lip, const ProbData<L, D> *const data) {
    AbsLoss<L, D>::Lip_initial_set(Lip,data);
    for (L row = 0; row < data->n; ++row)
        Lip[row]  += AbsLoss<L,D>::alpha;
}


#endif //LOSSFUN_ABSLOSS_H

