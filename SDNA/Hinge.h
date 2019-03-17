//
// Created by Lu Meng on 2018/5/22.
//

#ifndef SDNA_HINGE_H
#define SDNA_HINGE_H

#include "GenInfo.h"

template <typename L, typename D>
class Hinge_dual:public GenInfo<L,D> {
public:
    D mu = 1.;
    D lambda;

    Hinge_dual(const std::string &root_path_name, const std::string &data_file_name, const bool &apply_leverage);

    void setLowerUpper();
    void setv(const L &tau);
    void w_initial(std::vector<D> &w, const std::vector<D> &x)const;
    D PrimalDualGap(D &primal,D &dual, const std::vector<D> &x, const std::vector<D> &w)const;
    void Gsub(std::vector<D> &Gs, const L &tau, const std::vector<L> &S, const std::vector<D> &x,
              const std::vector<D> &w)const;
    void w_update(std::vector<D> &w, const L &tau, const std::vector<L> &S,const std::vector<D> &delta_x)const;

};

template<typename L, typename D>
Hinge_dual<L, D>::Hinge_dual(const std::string &root_path_name, const std::string &data_file_name, const bool &apply_leverage):
        GenInfo<L,D>(root_path_name, data_file_name, apply_leverage) {
    for (L row = 0; row < GenInfo<L,D>::n; ++row) {
        for (L col = GenInfo<L,D>::A_row_ptr[row]; col < GenInfo<L,D>::A_row_ptr[row+1]; ++col) {
            GenInfo<L,D>::A_value[col] *= GenInfo<L,D>::b[row];
        }
    }
    GenInfo<L,D>::setATbeA();
    L n = GenInfo<L,D>::n;
    L d = GenInfo<L,D>::d;
    Print("Calculate A..");
    Mtranspose(d,n,GenInfo<L,D>::AT_row_ptr,GenInfo<L,D>::AT_col_idx,GenInfo<L,D>::AT_value,
               GenInfo<L,D>::A_row_ptr,GenInfo<L,D>::A_col_idx,GenInfo<L,D>::A_value);
    GenInfo<L,D>::setATandAinfo();
    GenInfo<L,D>::setB(root_path_name,data_file_name,apply_leverage);
    lambda= param.lambda_multiple/d;
    GenInfo<L,D>::set_scale_add_L(1./(lambda*d*d),mu/d);
    GenInfo<L,D>::scalevalue = 1./(lambda*d*d);
    GenInfo<L,D>::lambda_min = mu/d;
    setLowerUpper();
}

template<typename L, typename D>
void Hinge_dual<L, D>::setLowerUpper() {
    L d = GenInfo<L,D>::d;
    GenInfo<L,D>::bound_idx = d;
    GenInfo<L,D>::lower_bound.clear();
    GenInfo<L,D>::upper_bound.clear();
    GenInfo<L,D>::lower_bound.resize(d,0);
    GenInfo<L,D>::upper_bound.resize(d,1);
}
// H = mu/d * I + 1/(lambda*d*d) *AT*A
template<typename L, typename D>
void Hinge_dual<L, D>::setv(const L &tau) {
    L d = GenInfo<L,D>::d;
    L n = GenInfo<L,D>::n;
    D beta_sum;
    GenInfo<L,D>::setv_prepare(beta_sum,tau);
    for (L row = 0; row < d; ++row) {
        GenInfo<L,D>::v[row] *= 1./(lambda*d*d);
        GenInfo<L,D>::v[row] += mu*beta_sum/(n*d);
    }
}

//w = 1/(lambda*d)Ax
template<typename L, typename D>
void Hinge_dual<L, D>::w_initial(std::vector<D> &w, const std::vector<D> &x)const {
    L d = GenInfo<L,D>::d;
    w = Ax(GenInfo<L,D>::n,GenInfo<L,D>::A_row_ptr,GenInfo<L,D>::A_col_idx,GenInfo<L,D>::A_value,x);
    w = AlphatimesV(d,1./(lambda*d),w);
}


/*!
    * Compute the gap between primal and dual.
    *
    * @param x : R^d
    * @param w : R^n
    * @param primal = 1/d * \sum phi_i(yi) + 0.5 * lambda *||w||^2
    * @param dual = -[1/d * (0.5 * mu * ||x||^2 - \sum xi) + 0.5 * lambda * ||w||^2]
    * @return gap = primal - dual
    * const guarantee this function will not update any members in class
    */
template<typename L, typename D>
D Hinge_dual<L, D>::PrimalDualGap(D &primal, D &dual, const std::vector<D> &x, const std::vector<D> &w) const {
    primal = 0;
    dual = 0;
    L d = GenInfo<L,D>::d;
    L n = GenInfo<L,D>::n;

    std::vector<D> y = Ax(d,GenInfo<L,D>::AT_row_ptr,GenInfo<L,D>::AT_col_idx,GenInfo<L,D>::AT_value,w);
    for (L row = 0; row < d; ++row) {
        if (y[row] < 1 - mu)
            primal += 1. - 0.5 * mu - y[row];
        else if (y[row] < 1)
            primal += 0.5 / mu * (1. - y[row]) * (1. - y[row]);
        else
            primal += 0.0;

        dual += 0.5 * mu * x[row] * x[row] - x[row];
    }
    D tmp = VtimesV(n,w,w);

    primal = primal/d + lambda * 0.5 * tmp;
    dual = -dual/d - lambda * 0.5 * tmp;

    return primal - dual;
}

//G = 1/d(mu*x - 1 + AT*w)
template<typename L, typename D>
void Hinge_dual<L, D>::Gsub(std::vector<D> &Gs, const L &tau, const std::vector<L> &S, const std::vector<D> &x,
                            const std::vector<D> &w) const {
    for (int row = 0; row < tau; ++row) {
        Gs[row] = mu * x[S[row]] - 1;
        for (L col = GenInfo<L, D>::AT_row_ptr[S[row]]; col < GenInfo<L, D>::AT_row_ptr[S[row] + 1]; ++col) {
            Gs[row] +=GenInfo<L,D>::AT_value[col] * w[GenInfo<L,D>::AT_col_idx[col]];
        }
        Gs[row] /= GenInfo<L,D>::d;
    }

}

//w = 1/(lambda*d) * Ax
template<typename L, typename D>
void Hinge_dual<L, D>::w_update(std::vector<D> &w, const L &tau, const std::vector<L> &S,
                                const std::vector<D> &delta_x) const {
    WplusAdelta(w,tau,S,GenInfo<L,D>::AT_row_ptr,GenInfo<L,D>::AT_col_idx,GenInfo<L,D>::AT_value,delta_x,
                1./(lambda*GenInfo<L,D>::d));

}


#endif //SDNA_HINGELOSS_H
