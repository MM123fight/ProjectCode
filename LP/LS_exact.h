//
// Created by Lu Meng on 2018/12/6.
//

#ifndef LP_LS_EXACT_H
#define LP_LS_EXACT_H

/*
 * x = (AJT*AJ+scale*I)^{-1}d = scale^{-1}(I-AJT*(AJ*AJT+scale*I)^{-1}AJ)d
 *   = 1/scale*(d - AJT*(AJ*AJT+scale*I)^{-1}AJd)
 */
template <typename L, typename D>
void LS_exact(std::vector<D> &x, D &time, const L& size, const L&n, const std::vector<L> &AJ_row_ptr, const std::vector<L> &AJ_col_idx,
              const std::vector<D> & AJ_value, const std::vector<D> &d, D scale = 1){
    D time_start;
    std::vector<D> HJ(size*size,0);
    std::vector<D> AJd(size,0);
    std::vector<D> sol(size,0);
    std::vector<D> AJTsol(n,0);
    gsl_matrix_view HJ_gsl = gsl_matrix_view_array(&HJ[0], size, size);
    gsl_vector_view AJd_gsl = gsl_vector_view_array(&AJd[0], size);
    gsl_vector_view T= gsl_vector_view_array(&sol[0], size);
    gsl_permutation * p = gsl_permutation_alloc(size);
    int s;

    time_start = clock();
    AtimesAT(HJ,size,n,AJ_row_ptr,AJ_col_idx,AJ_value);
    for (L row = 0; row < size; ++row)
        HJ[row*size+row] += scale;
    Ax(AJd,size,AJ_row_ptr,AJ_col_idx,AJ_value,d);
    gsl_linalg_LU_decomp(&HJ_gsl.matrix, p, &s);
    gsl_linalg_LU_solve(&HJ_gsl.matrix, p, &AJd_gsl.vector,&T.vector);
    ATx(AJTsol,size,n,AJ_row_ptr,AJ_col_idx,AJ_value,sol);
    for (L row = 0; row < n; ++row)
        x[row] = (d[row] - AJTsol[row])/scale;
    time += (D)(clock() - time_start)/CLOCKS_PER_SEC;
    gsl_permutation_free(p);
};
/* (AJT*AJ+scale*I)x = d;
*/


#endif //LP_LS_EXACT_H
