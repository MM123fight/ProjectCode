//
// Created by Lu Meng on 2018/5/10.
//

//Solve from primal
/* Primal problem: A is a n*d matrix,x: length d      Dual function value:(optimal: y = Ax - b )
 *
 * min{ 1/2||Ax - b||^2 }                            max -{ 1/2 y^2 + b^Ty}
 * s.t. x_[s] >=0,                                    s.t. [A^Ty]_[s] >= 0
 *                                                         [A^Ty]_[s~m] = 0
 *
 * KKT condition: A^T(Ax-b)_[s] >= 0, A^T(Ax-b)_[s~m] =0
 */

//First construct LP to LS problem
/* LP problem: A is a n*m matrix         can be transformed into a Least square problem with
 * A:n*d, b:R^n, x:R^d, c:R^d                    (size d)  (size d)  (size n)
 * min <c,x>                             A_LS =     0        I        A^T   (size d)  x_LS = x:R^d  b_LS  = c
 * Ax = b                                           A        0         0    (size n)         z:R^d          b
 * x >= 0                                          c^T       0       -b^T   (size 1)         y:R^n          0
 *                                       x >=0, z >=0(The first 2d are non-negative)
 *                                              (size d)  (size n)  (size 1)
 *                                       AT_LS =    0       A^T        c     (size d)
 *                                                  I        0         0     (size d)
 *                                                  A        0        -b     (size n)
 */

#ifndef SDNA_LS_H
#define SDNA_LS_H

#include "GenInfo.h"

template <typename L, typename D>
class LS_primal:public GenInfo<L,D>{
public:
    LS_primal(const std::string &root_path_name, const std::string &data_file_name, const bool &apply_leverage = false);
    virtual ~LS_primal() {}

    D PrimalDualGap(D &primal,D &dual, const std::vector<D> &x, const std::vector<D> &w)const;
    void Gsub(std::vector<D> &Gs,const L &tau, const std::vector<L> &S, const std::vector<D> &x, const std::vector<D> &w)const;
    void w_update(std::vector<D> &w, const L &tau, const std::vector<L> &S,const std::vector<D> &delta_x)const;

};

template <typename L, typename D>
class LPtoLS:public LS_primal<L,D> {
public:
    LPtoLS(const std::string &root_path_name, const std::string &data_file_name, const std::string &c_file_name);
    virtual ~LPtoLS(){}

};


template<typename L, typename D>
LS_primal<L, D>::LS_primal(const std::string &root_path_name, const std::string &data_file_name,
                           const bool &apply_leverage):GenInfo(root_path_name, data_file_name, apply_leverage) {
    L d = GenInfo<L,D>::d;
    L n = GenInfo<L,D>::n;
    Print("Calculate AT..");
    Mtranspose(n,d,GenInfo<L,D>::A_row_ptr,GenInfo<L,D>::A_col_idx,GenInfo<L,D>::A_value,GenInfo<L,D>::AT_row_ptr,GenInfo<L,D>::AT_col_idx,GenInfo<L,D>::AT_value);
    GenInfo<L,D>::setATandAinfo();
    GenInfo<L,D>::setB(root_path_name,data_file_name,apply_leverage);

}
//w = Ax - b
//primal = 1/2||w||^2,
template<typename L, typename D>
D LS_primal<L, D>::PrimalDualGap(D &primal, D &dual, const std::vector<D> &x, const std::vector<D> &w) const {
    L n = GenInfo<L,D>::n;
    primal = 0.5*l2norm_square(n,w);
    dual = 0;
    return primal;
}

//G = AT*(Ax-b) = AT*w
template<typename L, typename D>
void LS_primal<L, D>::Gsub(std::vector<D> &Gs, const L &tau, const std::vector<L> &S, const std::vector<D> &x,
                           const std::vector<D> &w) const {
    for (int row = 0; row < tau; ++row) {
        Gs[row] = 0;
        for (L col = GenInfo<L, D>::AT_row_ptr[S[row]]; col < GenInfo<L, D>::AT_row_ptr[S[row] + 1]; ++col) {
            Gs[row] +=GenInfo<L,D>::AT_value[col] * w[GenInfo<L,D>::AT_col_idx[col]];
        }
    }
}

// w = Ax-b = w+ A*delta_x;
template<typename L, typename D>
void LS_primal<L, D>::w_update(std::vector<D> &w, const L &tau, const std::vector<L> &S,
                                 const std::vector<D> &delta_x) const {
    WplusAdelta(w,tau,S,GenInfo<L,D>::AT_row_ptr,GenInfo<L,D>::AT_col_idx,GenInfo<L,D>::AT_value,delta_x,1.);
}



template<typename L, typename D>
LPtoLS<L, D>::LPtoLS(const std::string &root_path_name, const std::string &data_file_name, const std::string &c_file_name):
    LS_primal<L,D>(root_path_name, data_file_name) {
    L n = GenInfo<L,D>::n;
    L d = GenInfo<L,D>::d;
    GenInfo<L,D>::bound_idx = 2 * d;
    //Read vector c;
    std::vector<D> c;
    Print("Loading vector c...");
    std::string full_data_path = GenInfo<L,D>::root_path + "/" + c_file_name;
    VectorRead(full_data_path,d,c);
    // We first calculate AT_LS
    /*       (size d)  (size n)  (size 1)
      AT_LS =    0       A^T        c     (size d)
                 I        0         0     (size d)
                 A        0        -b     (size n)
    */
    Print("Calculate AT..");

    // A^T insert c
    L insert_index;
    for (L i = 0; i < d; ++i) {
        for (L j = GenInfo<L, D>::AT_row_ptr[i]; j < GenInfo<L, D>::AT_row_ptr[i + 1]; ++j)
            GenInfo<L, D>::AT_col_idx[j] += d;
    }
    for (L i = 0; i < d; ++i) {
        insert_index = GenInfo<L,D>::AT_row_ptr[i+1]+i;
        GenInfo<L,D>::AT_col_idx.insert(GenInfo<L,D>::AT_col_idx.begin() + insert_index, n+d);
        GenInfo<L,D>::AT_value.insert(GenInfo<L,D>::AT_value.begin() + insert_index, c[i]);
    }
    for (L i = 0; i < d; ++i)
        GenInfo<L, D>::AT_row_ptr[i + 1] += i + 1;

    // A^T c insert In
    std::vector<D> unitone(d,1);
    std::vector<L> index(d);
    std::vector<L> ptr(d);
    D last = 1+GenInfo<L,D>::AT_row_ptr.back();
    for (L i = 0; i < d; ++i) {
        index[i] = i;
        ptr[i] = i+last;
    }
    GenInfo<L,D>::AT_value.insert(GenInfo<L,D>::AT_value.end(), unitone.begin(),unitone.end());
    GenInfo<L,D>::AT_col_idx.insert(GenInfo<L,D>::AT_col_idx.end(), index.begin(),index.end());
    GenInfo<L,D>::AT_row_ptr.insert(GenInfo<L,D>::AT_row_ptr.end(), ptr.begin(),ptr.end());

    // A^T c In insert  A b
    L start;
    L end;

    for (L i = 0; i < n; ++i) {
        start = GenInfo<L,D>::A_row_ptr[i];
        end = GenInfo<L,D>::A_row_ptr[i+1];
        GenInfo<L,D>::AT_value.insert(GenInfo<L,D>::AT_value.end(),
                                          GenInfo<L,D>::A_value.begin()+start,
                                          GenInfo<L,D>::A_value.begin()+end);
        GenInfo<L,D>::AT_value.push_back(-GenInfo<L,D>::b[i]);
        GenInfo<L,D>::AT_col_idx.insert(GenInfo<L,D>::AT_col_idx.end(),
                                           GenInfo<L,D>::A_col_idx.begin()+start,
                                           GenInfo<L,D>::A_col_idx.begin()+end);
        GenInfo<L,D>::AT_col_idx.push_back(d+n);
        GenInfo<L,D>::AT_row_ptr.push_back(GenInfo<L,D>::AT_row_ptr.back()+end-start+1);
    }
    Print("Finish the Calculation of AT..");

    //Calculate A_LS
    GenInfo<L,D>::n = n+d+1;
    GenInfo<L,D>::d = n+2*d;
    d = GenInfo<L,D>::d;
    n = GenInfo<L,D>::n;
    Print("Calculate A..");
    Mtranspose(d,n, GenInfo<L,D>::AT_row_ptr,GenInfo<L,D>::AT_col_idx,GenInfo<L,D>::AT_value,
                          GenInfo<L,D>::A_row_ptr,GenInfo<L,D>::A_col_idx,GenInfo<L,D>::A_value);
    Print("Finish the Calculation of A..");

    //Caculate b_LS;
    GenInfo<L,D>::b.insert(GenInfo<L,D>::b.begin(),c.begin(),c.end());
    GenInfo<L,D>::b.push_back(0);

    GenInfo<L,D>::setATbeA();

}


#endif //SDNA_LEASTSQUARE_H
