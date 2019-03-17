//
// Created by Lu Meng on 2018/8/12.
//

#ifndef SDNA_L2LS_H
#define SDNA_L2LS_H

#include "GenInfo.h"
/*!
    * solve from primal
    * A:n*d, w : R^n, x : R^d
    * primal = 1/(2n) * ||Ax - b||^2 + lambda/2 *||x||^2
    *        = 1/(2n) * ||w||^2 + lambda/2 *||x||^2
    * dual = -[1/d * (1/2 * ||w||^2 - <w,b>) + lambda/2 * ||1/(lambda*d)A^Tw||^2]
    * w = b - Ax
    * H = 1/n*A^TA + lambda*I
    * G = 1/n*A^T(Ax-b) + lambda*x = lambda*x - 1/n*A^Tw
    */
template <typename L, typename D>
class L2LS_primal:public GenInfo<L,D> {
public:

    D lambda= param.lambda_multiple/GenInfo<L,D>::n;

    L2LS_primal(const std::string &root_path_name, const std::string &data_file_name, const bool &apply_leverage);
    virtual ~L2LS_primal(){};
    void setv(const L &tau);

    void w_initial(std::vector<D> &w, const std::vector<D> &x)const;
    D PrimalDualGap(D &primal,D &dual, const std::vector<D> &x, const std::vector<D> &w)const;
    void Gsub(std::vector<D> &Gs, const L &tau, const std::vector<L> &S,const std::vector<D> &x, const std::vector<D> &w)const;
    void w_update(std::vector<D> &w, const L &tau, const std::vector<L> &S,const std::vector<D> &delta_x)const;
};

/*!
    * solve from dual: First, you need to transform A to AT
    * A:n*d, w : R^n, x : R^d
    * primal = 1/(2d) * ||A^Tw - b||^2 + lambda/2 *||w||^2
    * dual = -[1/d * (1/2 * ||x||^2 - <x,b>) + lambda/2 * ||1/(lambda*d)Ax||^2]
    *      = -[1/d * (1/2 * ||x||^2 - <x,b>) + lambda/2 * ||w||^2]
    * w = 1/(lambda*d)Ax
    * H = 1/d*I +1/(lambda*d^2)A^TA
    * G = 1/d*(x-b+A^Tw)
    */
template <typename L, typename D>
class L2LS_dual:public GenInfo<L,D> {
public:

    D lambda;

    L2LS_dual(const std::string &root_path_name, const std::string &data_file_name, const bool &apply_leverage);
    virtual ~L2LS_dual() {}

    void setv(const L &tau);
    void w_initial(std::vector<D> &w, const std::vector<D> &x)const;
    D PrimalDualGap(D &primal,D &dual, const std::vector<D> &x, const std::vector<D> &w)const;
    void Gsub(std::vector<D> &Gs,const L &tau, const std::vector<L> &S, const std::vector<D> &x, const std::vector<D> &w)const;
    void w_update(std::vector<D> &w, const L &tau, const std::vector<L> &S,const std::vector<D> &delta_x)const;
};


template<typename L, typename D>
L2LS_primal<L, D>::L2LS_primal(const std::string &root_path_name, const std::string &data_file_name,const bool &apply_leverage):
        GenInfo<L,D>(root_path_name, data_file_name,apply_leverage) {
    L d = GenInfo<L,D>::d;
    L n = GenInfo<L,D>::n;
    Print("Calculate AT..");
    Mtranspose(n,d,GenInfo<L,D>::A_row_ptr,GenInfo<L,D>::A_col_idx,GenInfo<L,D>::A_value,GenInfo<L,D>::AT_row_ptr,GenInfo<L,D>::AT_col_idx,GenInfo<L,D>::AT_value);
    GenInfo<L,D>::setATandAinfo();
    GenInfo<L,D>::setB(root_path_name,data_file_name,apply_leverage);
    for (L row = 0; row < GenInfo<L,D>::d; ++row) {
        GenInfo<L, D>::ATL[row] *= 1./n;
        GenInfo<L, D>::ATL[row] += lambda;
    }
    for (L row = 0; row < GenInfo<L, D>::d; ++row) {
        GenInfo<L, D>::BTL[row] *= 1./n;
        GenInfo<L, D>::BTL[row] += lambda;
    }
    GenInfo<L,D>::scalevalue = 1./n;
}

template<typename L, typename D>
void L2LS_primal<L, D>::setv(const L &tau) {
    L d = GenInfo<L,D>::d;
    L n = GenInfo<L,D>::n;
    D beta_sum;
    GenInfo<L,D>::setv_prepare(beta_sum,tau);
    for (L row = 0; row < d; ++row) {
        GenInfo<L,D>::v[row] += beta_sum*lambda;
        GenInfo<L,D>::v[row] /= (D)n;
    }
}


//w = b - Ax
template<typename L, typename D>
void L2LS_primal<L, D>::w_initial(std::vector<D> &w, const std::vector<D> &x)const {
    L n = GenInfo<L,D>::n;
    w = Ax(n,GenInfo<L,D>::A_row_ptr,GenInfo<L,D>::A_col_idx,GenInfo<L,D>::A_value,x);
    w = VminusV(n,GenInfo<L,D>::b,w);
}



//primal= 1/(2n) * ||w||^2 + lambda/2 *||x||^2
//dual = -[1/n * (1/2 * ||w||^2 - <w,b>) + lambda/2 * ||1/(lambda*n)A^Tw||^2]
template<typename L, typename D>
D L2LS_primal<L, D>::PrimalDualGap(D &primal,D &dual, const std::vector<D> &x, const std::vector<D> &w) const {
    L d = GenInfo<L,D>::d;
    L n = GenInfo<L,D>::n;
    primal = 0;
    dual = 0;
    D tmp = 0.5/n*l2norm_square(n,w);
    primal = tmp + 0.5*lambda*(l2norm_square(d,x));
    dual = l2norm_square(d,Ax(d,GenInfo<L,D>::AT_row_ptr,GenInfo<L,D>::AT_col_idx,GenInfo<L,D>::AT_value,w));
    dual = -(tmp - 1./n*VtimesV(n,w,GenInfo<L,D>::b)+0.5/(lambda*n*n)*dual);
    return primal - dual;
}



//G = lambda*x - 1/n*A^Tw
template<typename L, typename D>
void L2LS_primal<L, D>::Gsub(std::vector<D> &Gs, const L &tau, const std::vector<L> &S, const std::vector<D> &x,
                                         const std::vector<D> &w) const {
    for (int row = 0; row < tau; ++row) {
        Gs[row] = 0;
        for (L col = GenInfo<L, D>::AT_row_ptr[S[row]]; col < GenInfo<L, D>::AT_row_ptr[S[row] + 1]; ++col) {
            Gs[row] +=GenInfo<L,D>::AT_value[col] * w[GenInfo<L,D>::AT_col_idx[col]];
        }
        Gs[row] = lambda*x[S[row]] - Gs[row]/GenInfo<L,D>::n;
    }
}

//w = b - Ax = w - A delta_x
template<typename L, typename D>
void L2LS_primal<L, D>::w_update(std::vector<D> &w, const L &tau, const std::vector<L> &S,
                                   const std::vector<D> &delta_x) const {
    WminusAdelta(w,tau,S,GenInfo<L,D>::AT_row_ptr,GenInfo<L,D>::AT_col_idx,GenInfo<L,D>::AT_value,delta_x,1.);
}


template<typename L, typename D>
L2LS_dual<L, D>::L2LS_dual(const std::string &root_path_name, const std::string &data_file_name, const bool &apply_leverage):
        GenInfo<L,D>(root_path_name, data_file_name, apply_leverage) {
    //Transform A to AT and AT to A
    GenInfo<L,D>::setATbeA();
    L d = GenInfo<L,D>::d;
    L n = GenInfo<L,D>::n;
    Print("Calculate A..");
    Mtranspose(d,n,GenInfo<L,D>::AT_row_ptr,GenInfo<L,D>::AT_col_idx,GenInfo<L,D>::AT_value,GenInfo<L,D>::A_row_ptr,GenInfo<L,D>::A_col_idx,GenInfo<L,D>::A_value);

    GenInfo<L,D>::setATandAinfo();
    GenInfo<L,D>::setB(root_path_name,data_file_name,apply_leverage);
    lambda= param.lambda_multiple/d;
    GenInfo<L,D>::set_scale_add_L(1./(lambda*d*d),1./d);
    GenInfo<L,D>::scalevalue = 1./(lambda*d*d);
}

// dual = -[1/d * (1/2 * ||x||^2 - <x,b>) + lambda/2 * ||1/(lambda*d)Ax||^2]
template<typename L, typename D>
void L2LS_dual<L, D>::setv(const L &tau) {
    L d = GenInfo<L,D>::d;
    L n = GenInfo<L,D>::n;
    D beta_sum;
    GenInfo<L,D>::setv_prepare(beta_sum,tau);
    for (L row = 0; row < d; ++row) {
        GenInfo<L,D>::v[row] *= 1./(lambda*d*d);
        GenInfo<L,D>::v[row] += beta_sum/(n*d);
    }
}


//w = 1/(lambda*d)Ax
template<typename L, typename D>
void L2LS_dual<L, D>::w_initial(std::vector<D> &w, const std::vector<D> &x)const {
    L n = GenInfo<L,D>::n;
    L d = GenInfo<L,D>::d;
    w = Ax(n,GenInfo<L,D>::A_row_ptr,GenInfo<L,D>::A_col_idx,GenInfo<L,D>::A_value,x);
    w = AlphatimesV(d,1./(lambda*d),w);
}

// primal = 1/(2d) * ||A^Tw - b||^2 + lambda/2 *||w||^2
// dual  = -[1/d * (1/2 * ||x||^2 - <x,b>) + lambda/2 * ||w||^2]
template<typename L, typename D>
D L2LS_dual<L, D>::PrimalDualGap(D &primal, D &dual, const std::vector<D> &x, const std::vector<D> &w) const {
    primal = 0;
    dual = 0;
    L d = GenInfo<L,D>::d;
    L n = GenInfo<L,D>::n;
    std::vector<D> tmp;
    D temp;
    temp = lambda/2*l2norm_square(n,w);
    tmp = Ax(d,GenInfo<L,D>::AT_row_ptr,GenInfo<L,D>::AT_col_idx,GenInfo<L,D>::AT_value,w);
    tmp = VminusV(d,tmp,GenInfo<L,D>::b);
    primal = 0.5/d*l2norm_square(d,tmp) + temp;
    dual = -(1./d*(0.5*l2norm_square(d,x) - VtimesV(d,x,GenInfo<L,D>::b)) + temp);
    return primal - dual;
}


//G = 1/d*(x-b+A^Tw)
template<typename L, typename D>
void L2LS_dual<L, D>::Gsub(std::vector<D> &Gs, const L &tau, const std::vector<L> &S, const std::vector<D> &x,
                                       const std::vector<D> &w) const {
    for (int row = 0; row < tau; ++row) {
        Gs[row] = x[S[row]] - GenInfo<L,D>::b[S[row]];
        for (L col = GenInfo<L, D>::AT_row_ptr[S[row]]; col < GenInfo<L, D>::AT_row_ptr[S[row] + 1]; ++col) {
            Gs[row] +=GenInfo<L,D>::AT_value[col] * w[GenInfo<L,D>::AT_col_idx[col]];
        }
        Gs[row] /= GenInfo<L,D>::d;
    }
}

// w = 1/(lambda*d)Ax
template<typename L, typename D>
void L2LS_dual<L, D>::w_update(std::vector<D> &w, const L &tau, const std::vector<L> &S,
                                 const std::vector<D> &delta_x) const {
    WplusAdelta(w,tau,S,GenInfo<L,D>::AT_row_ptr,GenInfo<L,D>::AT_col_idx,GenInfo<L,D>::AT_value,delta_x,
                1./(lambda*GenInfo<L,D>::d));
}



#endif //SDNA_L2LS_H
