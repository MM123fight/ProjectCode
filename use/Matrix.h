//
// Created by Lu Meng on 2018/10/23.
//

#ifndef USE_MATRIX_H
#define USE_MATRIX_H
#include "stdafx.h"
//Know the information of A and calculate the informmation of A^T
template <typename L, typename D>
void Mtranspose(std::vector<L> &MT_row_ptr, std::vector<L> &MT_col_idx,std::vector<D> &MT_value, const L &n, const L &d, const std::vector<L> &M_row_ptr, const std::vector<L> &M_col_idx,
                const std::vector<D> &M_value){

    MT_value.clear();
    MT_col_idx.clear();
    MT_row_ptr.clear();


    std::vector<std::vector<L>> tmp_idx(d);
    std::vector<std::vector<D>> tmp_value(d);
    for(L i = 0; i < n; ++i){
        for(L j = M_row_ptr[i]; j < M_row_ptr[i+1];++j){
            tmp_idx[M_col_idx[j]].push_back(i);
            tmp_value[M_col_idx[j]].push_back(M_value[j]);
        }
    }

    MT_row_ptr.resize(d+1,0);
    for(L i = 0; i < d; ++i){
        MT_value.insert(MT_value.end(), tmp_value[i].begin(), tmp_value[i].end());
        MT_col_idx.insert(MT_col_idx.end(), tmp_idx[i].begin(), tmp_idx[i].end());
        MT_row_ptr[i+1] = MT_col_idx.size();
    }

}

template <typename L, typename D>
void Mdiag(const L &n, const std::vector<L> &A_row_ptr,const std::vector<D> &A_value, std::vector<D> &AL) {
    AL.clear();
    AL.resize(n,0);
    for (L row = 0; row < n; ++row) {
        for (L col = A_row_ptr[row]; col < A_row_ptr[row+1]; ++col)
            AL[row] += A_value[col]*A_value[col];
    }
};

template <typename D>
inline void MtimesM_s2(std::vector<D> &A, const std::vector<D> &B) {
    std::vector<D> C(4);
    C[0] = A[0] * B[0] + A[1] * B[2];
    C[1] = A[0] * B[1] + A[1] * B[3];
    C[2] = A[2] * B[0] + A[3] * B[2];
    C[3] = A[2] * B[1] + A[3] * B[3];
    A[0] = C[0]; A[1] = C[1]; A[2] = C[2]; A[3] = C[3]; 

};

template <typename D>
inline void MtimesM_s2(std::vector<D>& C, const std::vector<D> &A, const std::vector<D> &B) {
    C[0] = A[0] * B[0] + A[1] * B[2];
    C[1] = A[0] * B[1] + A[1] * B[3];
    C[2] = A[2] * B[0] + A[3] * B[2];
    C[3] = A[2] * B[1] + A[3] * B[3];

};

template <typename L, typename D>
inline void MtimesM(std::vector<D>& C, const L& m, const L& n, const std::vector<D> &A, const std::vector<D> &B) {
    L idx1, idx2, idx3;
    for (L i = 0; i < m; ++i){
        for (L j= 0; j < n; ++j){
            idx1 = i*n + j;
            idx2 = i*n;
            idx3 = j*n;
            C[idx1] = 0;
            for (L row = 0; row < n; ++row)
                C[idx1] += A[idx2+row] * B[idx3+row];
        }
    }

};

//M = (alpha*I + beta*AA^T )
//AL:The diagonal matrix of M
template <typename L, typename D>
void Msub(std::vector<D> &Ms,const std::vector<L> &A_row_ptr, const std::vector<L> &A_col_idx,
          const std::vector<D> &A_value, const std::vector<D> &AL, const L &tau, const std::vector<L> &S,
                    const D &scalevalue = 1){
    Ms.clear();
    Ms.resize(tau*tau,0);

    D tmp = 0;
    L id1, id1_bound;
    L id2, id2_bound;
    for (L row = 0; row < tau; ++row) {
        for (L col = row+1; col < tau; ++col) {
            tmp = 0;
            id1 = A_row_ptr[S[row]], id1_bound = A_row_ptr[S[row] + 1];
            id2 = A_row_ptr[S[col]], id2_bound = A_row_ptr[S[col] + 1];
            while (id1 < id1_bound && id2 < id2_bound) {
                if (A_col_idx[id1] == A_col_idx[id2]) {
                    tmp += A_value[id1] * A_value[id2];
                    ++id1;
                    ++id2;
                } else if (A_col_idx[id1] < A_col_idx[id2]) {
                    ++id1;
                } else {
                    ++id2;
                }
            }
            tmp *= scalevalue;
            Ms[row * tau + col] = tmp;
            Ms[col * tau + row] = tmp;
        }
        Ms[row * tau + row] = AL[S[row]];
    }
};

template <typename L, typename D>
void AtimesAT(std::vector<D> &AAT, const L &m, const L &n, const std::vector<L> &A_row_ptr, const std::vector<L> &A_col_idx,
          const std::vector<D> &A_value){
    AAT.clear();
    AAT.resize(m*m,0);

    D tmp = 0;
    L id1, id1_bound;
    L id2, id2_bound;
    for (L row = 0; row < m; ++row) {
        for (L col = row; col < m; ++col) {
            tmp = 0;
            id1 = A_row_ptr[row], id1_bound = A_row_ptr[row + 1];
            id2 = A_row_ptr[col], id2_bound = A_row_ptr[col + 1];
            while (id1 < id1_bound && id2 < id2_bound) {
                if (A_col_idx[id1] == A_col_idx[id2]) {
                    tmp += A_value[id1] * A_value[id2];
                    ++id1;
                    ++id2;
                } else if (A_col_idx[id1] < A_col_idx[id2]) {
                    ++id1;
                } else {
                    ++id2;
                }
            }
            AAT[row * m + col] = tmp;
            AAT[col * m + row] = tmp;
        }
    }
    for (L row = 0; row < m; ++row){
        tmp = 0;
        for (L col = A_row_ptr[row]; col < A_row_ptr[row+1] ; ++col)
            tmp += A_value[col] * A_value[col];
        AAT[row*m+row] = tmp;
    }
};


//grad = ATw;
template <typename L, typename D>
inline void cord_ATw(std::vector<D> &grad, const L &tau, const std::vector<L> &S, const std::vector<L> &AT_row_ptr,
                        const std::vector<L> &AT_col_idx, const std::vector<D> & AT_value, const std::vector<D> w){
    L tmp;
    for (L row = 0; row < tau; ++row){
        tmp = S[row];
        grad[row] = 0;
        for (L col = AT_row_ptr[tmp]; col < AT_row_ptr[tmp+1]; ++col)
            grad[row] += AT_value[col] * w[AT_col_idx[col]];
    }
}

// w = w+alpha*A(delta_x)
template <typename L, typename D>
inline void WplusAdelta(std::vector<D> &w, const L &tau, const std::vector<L> &S, const std::vector<L> &AT_row_ptr,
                        const std::vector<L> &AT_col_idx, const std::vector<D> & AT_value, const std::vector<D> &delta,
                        const D &alpha){
    for (L row = 0; row < tau; ++row) {
        for (L col = AT_row_ptr[S[row]]; col < AT_row_ptr[S[row] + 1]; ++col)
            w[AT_col_idx[col]] += alpha*AT_value[col] * delta[row];
    }
};
// w = w-alpha*A(delta_x)
template <typename L, typename D>
inline void WminusAdelta(std::vector<D> &w, const L &tau, const std::vector<L> &S, const std::vector<L> &AT_row_ptr,
                        const std::vector<L> &AT_col_idx, const std::vector<D> & AT_value, const std::vector<D> &delta,
                         const D &alpha){
    for (L row = 0; row < tau; ++row) {
        for (L col = AT_row_ptr[S[row]]; col < AT_row_ptr[S[row] + 1]; ++col)
            w[AT_col_idx[col]] -= alpha*AT_value[col] * delta[row];
    }
};

// w = AT*x -b
template <typename L, typename D>
inline void ATxminusb(std::vector<D> &w, const L &m, const L& n, const std::vector<L> &A_row_ptr,
                      const std::vector<L> &A_col_idx, const std::vector<D> & A_value, const std::vector<D> &x,
                      const std::vector<D> &b){
    w.clear();
    w.resize(n,0);
    for (L row = 0; row < m; ++row) {
        for (L col = A_row_ptr[row]; col < A_row_ptr[row + 1]; ++col)
            w[A_col_idx[col]] += A_value[col] * x[row];
    }
    for (L row = 0; row < n; ++row)
        w[row] -= b[row];
};

template <typename L, typename D>
inline void ATx(std::vector<D> &w, const L &m, const L& n, const std::vector<L> &A_row_ptr,
                const std::vector<L> &A_col_idx, const std::vector<D> & A_value, const std::vector<D> &x){
    w.clear();
    w.resize(n,0);
    for (L row = 0; row < m; ++row) {
        for (L col = A_row_ptr[row]; col < A_row_ptr[row + 1]; ++col)
            w[A_col_idx[col]] += A_value[col] * x[row];
    }
};

template <typename L, typename D>
inline std::vector<D> Ax(const L &n, const std::vector<D> &H,const std::vector<D> &x){
    std::vector<D> y = std::vector<D>(n);
    for (L row = 0; row < n ; ++row) {
        for (L col = 0; col < n; ++col)
            y[row] += H[row*n+col] * x[col];
    }
    return y;
};

template <typename L, typename D>
inline std::vector<D> Ax(const L &n,const std::vector<L> &A_row_ptr, const std::vector<L> &A_col_idx,
                  const std::vector<D> & A_value,const std::vector<D> &x){
    std::vector<D> y = std::vector<D>(n);
    for (L row = 0; row < n ; ++row) {
        for (L col = A_row_ptr[row]; col < A_row_ptr[row+ 1]; ++col)
            y[row] += A_value[col] * x[A_col_idx[col]];
    }
    return y;
};

template <typename L, typename D>
inline void Ax(std::vector<D> &y, const L &n,const std::vector<L> &A_row_ptr, const std::vector<L> &A_col_idx,
                         const std::vector<D> & A_value,const std::vector<D> &x){
    y.clear();
    y.resize(n,0);
    for (L row = 0; row < n ; ++row) {
        for (L col = A_row_ptr[row]; col < A_row_ptr[row+ 1]; ++col)
            y[row] += A_value[col] * x[A_col_idx[col]];
    }
};

template <typename L, typename D>
inline void Axminusb(std::vector<D> &y, const L &n,const std::vector<L> &A_row_ptr, const std::vector<L> &A_col_idx,
               const std::vector<D> & A_value,const std::vector<D> &x, const std::vector<D> &b){
    y.clear();
    y.resize(n,0);
    for (L row = 0; row < n ; ++row) {
        for (L col = A_row_ptr[row]; col < A_row_ptr[row+ 1]; ++col)
            y[row] += A_value[col] * x[A_col_idx[col]];
        y[row] -= b[row];
    }
};

template <typename L, typename D>
inline void Ax(std::vector<D> &y, const L &n,const std::vector<D> &H, const std::vector<D> &x){
    y.clear();
    y.resize(n,0);
    for (L row = 0; row < n ; ++row) {
        for (L col = 0; col < n; ++col)
            y[row] += H[row*n+col] * x[col];
    }
};



template <typename L, typename D>
inline std::vector<D> VplusV(const L &n, const std::vector<D> &x , const std::vector<D> &y ){
    std::vector<D> z = std::vector<D>(n);
    for (L row = 0; row < n; ++row)
        z[row] = x[row] + y[row];
    return z;
};

template <typename L, typename D>
inline std::vector<D> VminusV(const L &n, const std::vector<D> &x , const std::vector<D> &y ){
    std::vector<D> z = std::vector<D>(n);
    for (L row = 0; row < n; ++row)
        z[row] = x[row] - y[row];
    return z;
};

template <typename L, typename D>
inline D VtimesV(const L &n, const std::vector<D> &x, const std::vector<D> &y){
    D tmp = 0;
    for (L row = 0;  row< n; ++row)
        tmp += x[row]*y[row];
    return tmp;
};

template <typename L, typename D>
inline std::vector<D> AlphatimesV(const L &n, const D &alpha, const std::vector<D> &y){
    std::vector<D> x(n);
    for (L row = 0;  row< n; ++row)
        x[row] = alpha*y[row];
    return x;
};

template <typename L,typename D>
inline D l2norm_square(const std::vector<D> &x){
    D tmp = 0;
    for (L row = 0;  row< x.size(); ++row)
        tmp += x[row]*x[row];
    return tmp;
};

template <typename L, typename D>
inline D l2norm_square(const L &n, const std::vector<D> &x){
    D tmp = 0;
    for (L row = 0;  row< n; ++row)
        tmp += x[row]*x[row];
    return tmp;
};

template <typename D>
void transform(std::vector<D> & x, std::vector<D> &y){
    std::vector<D> tmp = x;
    x = y;
    y = tmp;
};
template <typename D>
void transform(D &x, D &y){
    D tmp = x;
    x = y;
    y = tmp;
};

//A^TA <= B^TB
template <typename L, typename D>
void LeverageA(L &Bn, std::vector<L> &B_row_ptr, std::vector<L> &B_col_idx, std::vector<D> &B_value, const L &n,
               const L &d, const std::vector<L> &A_row_ptr, const std::vector<L> &A_col_idx, const std::vector<D> &A_value,
               const std::vector<D> &Leverage, const D &epsilon, const D &c){
    B_row_ptr.clear();
    B_col_idx.clear();
    B_value.clear();
    B_row_ptr.push_back(0);
    Bn = 0;

    D p;
    D tmp;
    L start;
    L end;
    D beta = pow(epsilon,-2)*c*log(d);
    Print(beta);
    std::default_random_engine e;
    for (L row = 0; row < n; ++row) {
        p = std::min(1.,beta*Leverage[row]);
        std::bernoulli_distribution b(p);
        if(b(e)){
            ++Bn;
            start = A_row_ptr[row];
            end = A_row_ptr[row+1];
            B_row_ptr.push_back(B_row_ptr.back()+ end - start);
            B_col_idx.insert(B_col_idx.end(),A_col_idx.begin()+start,A_col_idx.begin()+end);
            tmp = sqrt(p*(1-epsilon));
            for (L col = start; col < end; ++col)
                B_value.push_back(A_value[col]/tmp);
        }
    }

};

template <typename L, typename D>
void leverage_sample(std::vector<L> & Sample, const L &n, const L &d, const std::vector<D> &Leverage, const D &epsilon, const D &c){
    D p;
    D beta = pow(epsilon,-2)*c*log(d);
    Sample.clear();
    std::default_random_engine e;
    for (L row = 0; row < n; ++row) {
        p = std::min(1.,beta*Leverage[row]);
        std::bernoulli_distribution b(p);
        if(b(e)){
            Sample.push_back(row);
        }
    }

};

template <typename L, typename D>
D findMax(std::vector<D> &v) {
    L n = v.size();
    D max = v[0];
    for (L row = 0; row < n; ++row) {
        if( max < v[row])
            max = v[row];
    }
    return max;
}

template < typename L, typename D >
void subMat(std::vector<L> &AJ_row_ptr,std::vector<L> &AJ_col_idx, std::vector<D> &AJ_value,
            const std::vector<L> &A_row_ptr, const std::vector<L> &A_col_idx, const std::vector<D> &A_value,
            const L &sizeJ, const std::vector<L> &J){
    AJ_row_ptr.clear();
    AJ_col_idx.clear();
    AJ_value.clear();
    L tmp;
    L start;
    L end;
    AJ_row_ptr.resize(sizeJ+1,0);
    for (L row = 0; row < sizeJ; ++row) {
        tmp = J[row];
        start = A_row_ptr[tmp];
        end = A_row_ptr[tmp+1];
        AJ_row_ptr[row+1] = AJ_row_ptr[row] + end - start;
        AJ_col_idx.insert(AJ_col_idx.end(), A_col_idx.begin()+start, A_col_idx.begin()+end);
        AJ_value.insert(AJ_value.end(), A_value.begin()+start, A_value.begin()+end);
    }
}

template <typename L, typename D>
void subVec(std::vector<D> &vJ, const std::vector<D> &v, const L &sizeJ, const std::vector<L> &J){
    vJ.resize(sizeJ,0);
    for (L row = 0; row < sizeJ; ++row)
        vJ[row] = v[J[row]];
}

template <typename L, typename D>
void Diag_ATA(std::vector<D> &Diag,const L&m, const L&n, const std::vector<L> &A_row_ptr, const std::vector<L> &A_col_idx, const std::vector<D> &A_value){
    Diag.clear();
    Diag.resize(n,0);
    for (L row = 0; row < m; ++row) {
        for (L col = A_row_ptr[row]; col < A_row_ptr[row+1]; ++col)
            Diag[A_col_idx[col]] += A_value[col]*A_value[col];
    }
};

void gsl_matrix_inv(gsl_matrix *A, gsl_matrix*inverse)
{
    size_t n = A->size1;
    gsl_matrix *tmpA = gsl_matrix_alloc(n, n);
    gsl_matrix_memcpy(tmpA, A);
    gsl_permutation *p = gsl_permutation_alloc(n);
    int sign = 0;
    gsl_linalg_LU_decomp(tmpA, p, &sign);
    //inverse = gsl_matrix_alloc(n, n);
    gsl_linalg_LU_invert(tmpA, p, inverse);
    gsl_permutation_free(p);
    gsl_matrix_free(tmpA);
}

#endif //USE_MATRIX_H
