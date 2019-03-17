//
// Created by Lu Meng on 2018/12/6.
//

#ifndef LP_CG_H
#define LP_CG_H
//Ax = b
template <typename L, typename D>
void CG(std::vector<D> &x, const L& size, const std::vector<D> &A, const std::vector<D> &b, const D& tol){
    std::vector<D> r(size,0);
    Ax(r,size,A,x);
    for (L row = 0; row < size; ++row)
        r[row] = b[row] - r[row];
    std::vector<D> p = r;
    std::vector<D> Ap;
    D r_square;
    D r_square_before;
    D pAp;
    D alpha;
    D beta;
    r_square = 0;
    for (L row = 0; row < size; ++row)
        r_square += r[row] * r[row];
    while(r_square > tol) {
        r_square_before = r_square;
        Ax(Ap, size, A, p);
        pAp = 0;
        for (L row = 0; row < size; ++row)
            pAp += p[row] * Ap[row];
        alpha = r_square_before/pAp;
        for (L row = 0; row < size; ++row) {
            x[row] += alpha * p[row];
            r[row] -= alpha * Ap[row];
        }
        r_square = 0;
        for (L row = 0; row < size; ++row)
            r_square += r[row] * r[row];
        beta = r_square/r_square_before;
        for (L row = 0; row < size; ++row)
            p[row] = r[row] + beta*p[row];
    }

};
/* (AJT*AJ+scale*I)x = d;
*/
template <typename L, typename D>
void CG(std::vector<D> &x, D &time, const L& size, const L&n, const std::vector<L> &AJ_row_ptr, const std::vector<L> &AJ_col_idx,
        const std::vector<D> & AJ_value, const std::vector<D> &d, const D& tol, D scale = 1){
    D time_start;
    std::vector<D> r(n,0);
    std::vector<D> AJx(size,0);
    std::vector<D> AJp(size,0);
    //(AJT*AJ+I)p
    std::vector<D> Hp;
    D r_square;
    D r_square_before;
    D pHp;
    D alpha;
    D beta;

    time_start = clock();
    Ax(AJx,size,AJ_row_ptr,AJ_col_idx,AJ_value,x);
    ATx(r,size,n,AJ_row_ptr,AJ_col_idx,AJ_value,AJx);
    for (L row = 0; row < n; ++row) {
        r[row] += scale * x[row];
        r[row] = d[row] - r[row];
    }
    std::vector<D> p = r;
    r_square = l2norm_square(n,r);
    while(r_square > tol) {
        r_square_before = r_square;
        Ax(AJp,size,AJ_row_ptr,AJ_col_idx,AJ_value,p);
        ATx(Hp,size,n,AJ_row_ptr,AJ_col_idx,AJ_value,AJp);
        for (L row = 0; row < n; ++row)
            Hp[row] += scale * p[row];
        pHp = 0;
        for (L row = 0; row < n; ++row)
            pHp += p[row] * Hp[row];
        alpha = r_square_before/pHp;
        for (L row = 0; row < n; ++row) {
            x[row] += alpha * p[row];
            r[row] -= alpha * Hp[row];
        }
        r_square = l2norm_square(n,r);
        beta = r_square/r_square_before;
        for (L row = 0; row < n; ++row)
            p[row] = r[row] + beta*p[row];
    }
    time += (D)(clock() - time_start)/CLOCKS_PER_SEC;

};
/*
 * x = (AJT*AJ+scale*I)^{-1}d = scale^{-1}(I-AJT*(AJ*AJT+scale*I)^{-1}AJ)d
 *   = 1/scale*(d - AJT*(AJ*AJT+scale*I)^{-1}AJd)
 *(AJ*AJT+scale*I)z = AJd
 * x = 1/scale(d - AJT*z)
 * (AJT*AJ+scale*I)x -d = 1/scale(AJT*AJ+scale*I)(d - AJT*z) -d  = 1/scale(AJT*AJ*d -AJT*(AJ*AJT+scale*I)*z)
 *                = 1/scale*AJT*(AJ*d - (AJ*AJT+scale*I)*z) = 1/scale*AJT*r;
 */

template <typename L, typename D>
void CG_inv(std::vector<D> &x, D &time, const L& size, const L&n, const std::vector<L> &AJ_row_ptr, const std::vector<L> &AJ_col_idx,
        const std::vector<D> & AJ_value, const std::vector<D> &d, const D& tol, D scale = 1) {
    D time_start;
    D inv_tol = pow(scale, 2) * tol;
    std::vector <D> AJd(size, 0);
    std::vector <D> z(size, 0);
    std::vector <D> r(size, 0);
    std::vector <D> AJTz(n, 0);
    std::vector <D> AJTp(n, 0);
    std::vector <D> Hp;
    D r_square;
    D r_square_before;
    D res_square;
    D pHp;
    D alpha;
    D beta;

    time_start = clock();
    Ax(AJd, size, AJ_row_ptr, AJ_col_idx, AJ_value, d);
    ATx(AJTz, size, n, AJ_row_ptr, AJ_col_idx, AJ_value, z);
    Ax(r, size, AJ_row_ptr, AJ_col_idx, AJ_value, AJTz);
    for (L row = 0; row < size; ++row) {
        r[row] += scale * z[row];
        r[row] = AJd[row] - r[row];
    }
    std::vector <D> p = r;
    ATx(AJTp, size, n, AJ_row_ptr, AJ_col_idx, AJ_value, p);
    std::vector <D> res = AJTp;
    r_square = l2norm_square(size, r);
    res_square = l2norm_square(n, res);

    while(res_square > inv_tol) {
        r_square_before = r_square;
        Ax(Hp, size, AJ_row_ptr, AJ_col_idx, AJ_value, AJTp);
        for (L row = 0; row < size; ++row)
            Hp[row] += scale * p[row];
        pHp = 0;
        for (L row = 0; row < size; ++row)
            pHp += p[row] * Hp[row];
        alpha = r_square_before / pHp;
        for (L row = 0; row < size; ++row) {
            z[row] += alpha * p[row];
            r[row] -= alpha * Hp[row];
        }
        r_square = l2norm_square(size, r);
        beta = r_square / r_square_before;
        for (L row = 0; row < size; ++row)
            p[row] = r[row] + beta * p[row];
        ATx(res, size, n, AJ_row_ptr, AJ_col_idx, AJ_value, r);
        for (L row = 0; row < n; ++row)
            AJTp[row] = res[row] + beta * AJTp[row];
        res_square = l2norm_square(n, res);
    }
    ATx(AJTz,size,n,AJ_row_ptr,AJ_col_idx,AJ_value,z);
    for (L row = 0; row < n; ++row)
        x[row] = (d[row] - AJTz[row])/scale;
    time += (D)(clock() - time_start)/CLOCKS_PER_SEC;

};
#endif //LP_CG_H
