//
// Created by Lu Meng on 2019/2/20.
//

#ifndef LOSSFUN_ABSLOSS_H
#define LOSSFUN_ABSLOSS_H

#include "../Problem/ProblemHeader.h"
#include "../use/useheader.h"

// w = Ax - b

template<typename L, typename D>
class AbsLoss{
protected:
public:

    D alpha = 0;
    std::vector<D> d = std::vector<D>();
    std::vector<D> b = std::vector<D>();


    AbsLoss(){}
    AbsLoss(const std::vector<D>& b_value): b(b_value){}
    AbsLoss(const D& alpha_value,const std::vector<D>& b_value, const std::vector<D>& d_value):
            alpha(alpha_value),b(b_value),d(d_value){}

    virtual ~AbsLoss() {}


    //w = Ax -b
    void w_initial(std::vector<D> &w, const ProbData<L, D>* const data, const std::vector<D> &x);
    void w_initial(std::vector<D> &w, const ProbData<L, D>* const data, const std::vector<D> &x,
                   const std::vector<D>& b);
    //w_dual  = AT*lambda+c
    void w_dual_initial(std::vector<D> &w_dual, const ProbData<L, D>* const data, const std::vector<D> &lambda,
                   const std::vector<D>& c);
    void w_update(std::vector<D> &w, const ProbData<L, D> *const data, const std::vector<D>& delta,
                  const L &blocksize, const std::vector<L> &cord);

    void set_b(const std::vector<D> &b);
    void set_d(const std::vector<D> &d);


    virtual void cord_grad_update(std::vector<D>& cord_grad, const ProbData<L, D>* const data, const std::vector<D>& w,
                                  const std::vector<D>& x, const L& blocksize, const std::vector<L>& cord){}
    virtual void grad_compute(std::vector<D>& grad, const ProbData<L, D>* const data,const std::vector<D>& w,
                              const std::vector<D>& x){}
    virtual void cord_grad_update(D& cord_grad, const ProbData<L, D>* const data, const D& coeff_1, const D& coeff_2,
                                  const std::vector<D>& w_y_bar, const std::vector<D>& y_bar, const std::vector<D>& w_z_bar,
                                  const std::vector<D>& z_bar, const L& cord){}
    virtual void grad_compute(std::vector<D> &grad, const ProbData<L, D>* const data, const D &tau,
                      const std::vector<D> &M, const std::vector<D> &w_y_bar, const std::vector<D> &y_bar,
                      const std::vector<D> &w_z_bar, const std::vector<D> &z_bar){}
    virtual void cord_grad_sum_update(std::vector<D>& cord_grad_sum, const ProbData<L, D>* const data, const std::vector<D>& w_y,
                                      const std::vector<D>& y, const std::vector<D>& w_z, const std::vector<D>& z,
                                      const D& theta_y, const L& blocksize, const std::vector<L>& cord){}
    virtual void grad_sum_compute(std::vector<D> &grad_sum,const ProbData<L, D>* const data, const std::vector<D> &w_y,
                                  const std::vector<D> &y, const std::vector<D> &w_z, const std::vector<D> &z, const D &theta_y){}
    virtual void fun_value_compute(D& fun_value, const ProbData<L, D>* const data, const std::vector<D>& w, const std::vector<D>& x){}
    //0.5*||Ax-b||~ESO(v)
    void v_initial_set(std::vector<D>& v, const ProbData<L, D>* const data, const L& blocksize);
    //0.5*||Ax-b||+g(x)~ESO(v)
    virtual void v_set(std::vector<D>& v, const ProbData<L, D>* const data, const L& blocksize){}
    void Lip_initial_set(std::vector<D>& Lip, const ProbData<L, D>* const data);
    virtual void Lip_set(std::vector<D>& Lip, const ProbData<L, D>* const data){}


};

template<typename L, typename D>
void AbsLoss<L, D>::w_initial(std::vector<D> &w,const ProbData<L, D>* const data,const std::vector<D> &x){
    for (L row = 0; row < data->m; ++row) {
        w[row] = 0;
        for (L col = data->A_row_ptr[row]; col < data->A_row_ptr[row+1]; ++col)
            w[row] += data->A_value[col] * x[data->A_col_idx[col]];
    }
}

template<typename L, typename D>
void AbsLoss<L, D>::w_initial(std::vector<D> &w,const ProbData<L, D>* const data,const std::vector<D> &x,
                             const std::vector<D>& b){
    for (L row = 0; row < data->m; ++row) {
        w[row] = -b[row];
        for (L col = data->A_row_ptr[row]; col < data->A_row_ptr[row+1]; ++col)
            w[row] += data->A_value[col] * x[data->A_col_idx[col]];
    }
}

template<typename L, typename D>
void AbsLoss<L, D>::w_dual_initial(std::vector<D> &w_dual, const ProbData<L, D>* const data, const std::vector<D> &lambda,
                      const std::vector<D>& c){
    for (L row = 0; row < data->n; ++row) {
        w_dual[row] = c[row];
        for (L col = data->AT_row_ptr[row]; col < data->AT_row_ptr[row+1]; ++col)
            w_dual[row] += data->AT_value[col] * lambda[data->AT_col_idx[col]];
    }
}
template<typename L, typename D>
void AbsLoss<L, D>::w_update(std::vector<D> &w, const ProbData<L, D> *const data, const std::vector<D>& delta, const L &blocksize,
                             const std::vector<L> &cord){
    for(L row = 0; row < blocksize; ++row) {
        for (L col = data->AT_row_ptr[cord[row]]; col< data->AT_row_ptr[cord[row] + 1]; ++col)
            w[data->AT_col_idx[col]] += delta[row] * data->AT_value[col];
    }
}

template<typename L, typename D>
void AbsLoss<L, D>::v_initial_set(std::vector<D> &v, const ProbData<L, D> *const data, const L& blocksize) {
    L n = data->n;
    L m = data->m;
    D eta = (D)(blocksize-1)/(D)std::max<L>(1,n - 1);
    std::vector<D> beta(m);
    D beta_sum = 0;
    for (L row = 0; row < m; ++row) {
        beta[row] = 1 + eta * (data->A_row_ptr[row+1] - data->A_row_ptr[row] - 1);
        beta_sum += beta[row];
    }
    v.clear();
    v.resize(n,0);
    for (L row = 0; row < n; ++row) {
        for (L col = data->AT_row_ptr[row]; col < data->AT_row_ptr[row + 1]; ++col)
            v[row] += beta[data->AT_col_idx[col]] * data->AT_value[col] * data->AT_value[col];
    }
}
template<typename L, typename D>
void AbsLoss<L, D>::Lip_initial_set(std::vector<D> &Lip, const ProbData<L, D> *const data) {
    for (L row = 0; row < data->n; ++row) {
        Lip[row] = 0;
        for (L col = data->AT_row_ptr[row]; col < data->AT_row_ptr[row + 1]; ++col)
            Lip[row] += data->AT_value[col] * data->AT_value[col];
    }

}

template<typename L, typename D>
void AbsLoss<L, D>::set_b(const std::vector<D> &b) {
        this->b = b;
}
template<typename L, typename D>
void AbsLoss<L, D>::set_d(const std::vector<D> &d) {
        this->d = d;
}



#endif //LOSSFUN_ABSLOSS_H

