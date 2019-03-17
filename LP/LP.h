//
// Created by Lu Meng on 2018/5/10.
//

#ifndef LP_LP_H
#define LP_LP_H

#include "../use/useheader.h"
#include "../SDNA/ProbParam.h"
#include "../SDNA/Print.h"
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_cblas.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include "LPuse.h"
#include "Param.h"
#include "NU_ACDM.h"
#include "SSN.h"


template <typename L, typename D>
class LP{
public:

    L mi,me,m;

    L nb,nf,n;

    std::vector<D> c = std::vector<D>();

    std::vector<L> Ae_row_ptr = std::vector<L>();
    std::vector<L> Ae_col_idx = std::vector<L>();
    std::vector<D> Ae_value = std::vector<D>();
    std::vector<D> be = std::vector<D>();

    std::vector<L> Ai_row_ptr = std::vector<L>();
    std::vector<L> Ai_col_idx = std::vector<L>();
    std::vector<D> Ai_value = std::vector<D>();
    std::vector<D> bi = std::vector<D>();

    L n_bar;
    L nb_bar;
    L m_bar;
    std::vector<L> A_row_ptr = std::vector<L>();
    std::vector<L> A_col_idx = std::vector<L>();
    std::vector<D> A_value = std::vector<D>();

    std::vector<L> AT_row_ptr = std::vector<L>();
    std::vector<L> AT_col_idx = std::vector<L>();
    std::vector<D> AT_value = std::vector<D>();

    std::vector<D> b = std::vector<D>();

    void setA1();
    void setA2();
    void setA3();
    void setAi();
    void setb1();
    void setb2();
    void setbi();
    std::vector<D> ADMM1(const D &tol,const std::string &subtype = "NU_ACDM");
    std::vector<D> ADMM2(const D &tol,const std::string &subtype = "NU_ACDM");
    std::vector<D> ADMM(const D &tol);
    //std::vector<D> SSNPPA(const D &tol);
    std::string root_path;
    std::ofstream logFile;

    LP(const std::string &root_path_name, const std::string &data_file_name);

    virtual ~LP() {}
};

template<typename L, typename D>
LP<L, D>::LP(const std::string &root_path_name, const std::string &data_file_name){
    root_path = root_path_name;
    std::string data_path = root_path_name + "/"+ data_file_name;
    Print("Loading Data...");
    readMeta(data_path+"/meta", mi, me, nb, nf);
    m = mi+me;
    n = nb+nf;
    readMat(data_path+"/A", mi, n,Ai_row_ptr, Ai_col_idx, Ai_value);
    readMat(data_path+"/Aeq", me, n,Ae_row_ptr, Ae_col_idx, Ae_value);
    VectorRead(data_path+"/b",mi,bi);
    VectorRead(data_path+"/beq",me,be);
    VectorRead(data_path+"/c",n,c);
    Print("Finish loading data from:", data_path);

}

template<typename L, typename D>
void LP<L, D>::setA1() {
    n_bar = n + mi;
    A_col_idx = Ae_col_idx;
    A_value = Ae_value;
    A_row_ptr = Ae_row_ptr;
    L tmp = A_row_ptr.back();
    //Print("tmp",tmp);
    //Print("Ai_row_ptr",Ai_row_ptr);
    for (L row = 0; row < mi; ++row) {
        A_col_idx.insert(A_col_idx.end(),Ai_col_idx.begin()+Ai_row_ptr[row],Ai_col_idx.begin()+Ai_row_ptr[row+1]);
        A_col_idx.push_back(n+row);
        A_value.insert(A_value.end(),Ai_value.begin()+Ai_row_ptr[row],Ai_value.begin()+Ai_row_ptr[row+1]);
        A_value.push_back(1);
        A_row_ptr.push_back(tmp+Ai_row_ptr[row+1]+row+1);
    }

}

template<typename L, typename D>
void LP<L, D>::setA2() {
    m_bar = 2*me + mi;
    A_col_idx = Ai_col_idx;
    A_value = Ai_value;
    A_row_ptr = Ai_row_ptr;
    A_col_idx.insert(A_col_idx.end(),Ae_col_idx.begin(),Ae_col_idx.end());
    A_value.insert(A_value.end(),Ae_value.begin(),Ae_value.end());
    A_col_idx.insert(A_col_idx.end(),Ae_col_idx.begin(),Ae_col_idx.end());
    for (L row = 0; row < Ae_row_ptr.back(); ++row)
        A_value.push_back(-Ae_value[row]);
    L tmp;
    for (int i = 0; i < 2; ++i) {
        tmp = A_row_ptr.back();
        for (L row = 0; row < me; ++row)
            A_row_ptr.push_back(tmp + Ae_row_ptr[row + 1]);
    }
}

template<typename L, typename D>
void LP<L, D>::setA3() {
    n_bar = n + mi;
    nb_bar = nb+mi;
    A_col_idx = Ae_col_idx;
    for (L row = 0; row < Ae_row_ptr.back(); ++row)
        A_col_idx[row] += mi;
    A_value = Ae_value;
    A_row_ptr = Ae_row_ptr;
    L tmp = A_row_ptr.back();
    std::vector<L> idx_tmp = Ai_col_idx;
    for (L row = 0; row < idx_tmp.size(); ++row)
        idx_tmp[row] +=mi;
    for (L row = 0; row < mi; ++row) {
        A_col_idx.push_back(row);
        A_col_idx.insert(A_col_idx.end(),idx_tmp.begin()+Ai_row_ptr[row],idx_tmp.begin()+Ai_row_ptr[row+1]);
        A_value.push_back(1);
        A_value.insert(A_value.end(),Ai_value.begin()+Ai_row_ptr[row],Ai_value.begin()+Ai_row_ptr[row+1]);
        A_row_ptr.push_back(tmp+Ai_row_ptr[row+1]+row+1);
    }

    //Print(m,n_bar,A_row_ptr,A_col_idx,A_value);

}

template<typename L, typename D>
void LP<L, D>::setAi() {
    std::vector<D> negones(nb,-1);
    Ai_value.insert(Ai_value.end(),negones.begin(),negones.end());
    L tmp = Ai_row_ptr.back()+1;
    for (L row = 0; row < nb; ++row) {
        Ai_row_ptr.push_back(tmp+row);
        Ai_col_idx.push_back(row);
    }
    mi = mi + nb;
    m = me +mi;
    A_row_ptr = Ae_row_ptr;
    A_col_idx = Ae_col_idx;
    A_value = Ae_value;
    A_value.insert(A_value.end(),Ai_value.begin(),Ai_value.end());
    A_col_idx.insert(A_col_idx.end(),Ai_col_idx.begin(),Ai_col_idx.end());
    tmp = A_row_ptr.back();
    for (L row = 0; row < mi; ++row)
        A_row_ptr.push_back(tmp+Ai_row_ptr[row+1]);
}

template<typename L, typename D>
void LP<L, D>::setb1() {
    b = be;
    b.insert(b.end(),bi.begin(),bi.end());
}

template<typename L, typename D>
void LP<L, D>::setb2() {
    b = bi;
    b.insert(b.end(),be.begin(),be.end());
    for (L row = 0; row < me; ++row)
        b.push_back(-be[row]);
}

template<typename L, typename D>
void LP<L, D>::setbi() {
    std::vector<D> zeros(nb,0);
    bi.insert(bi.end(),zeros.begin(),zeros.end());
}

template<typename L, typename D>
std::vector<D> LP<L, D>::ADMM1(const D &tol,const std::string &subtype) {
    D rho = lp_param.rho;
    setA1();
    setb1();
    D time = 0;
    D Matrix_time;
    D time_start;
    std::vector<D> x(n_bar,0);
    std::vector<D> y(n_bar,0);
    std::vector<D> lambda_E(m,0);
    std::vector<D> lambda_x(n_bar,0);
    std::vector<D> Atimesx(m,0);
    std::vector<D> AATplusI(m*m,0);
    std::vector<D> Inv(m*m,0);
    gsl_matrix_view H = gsl_matrix_view_array(&AATplusI[0],m,m);
    gsl_matrix_view Hinv = gsl_matrix_view_array(&Inv[0],m,m);

    std::vector<D> we(m,0);
    std::vector<D> ve(n_bar,0);
    std::vector<D> d(m,0);
    std::vector<D> xminusy(n_bar,0);
    std::vector<D> Aeveminuswe(m,0);
    std::vector<D> g(m,0);
    std::vector<D> ge(n_bar,0);
    std::vector<D> dual_res(n_bar,0);
    std::vector<D> dual_res_E(n_bar,0);
    std::vector<D> dual_res_x(n_bar,0);
    std::vector<D> ATb(n_bar,0);
    ATx(ATb,m,n_bar,A_row_ptr,A_col_idx,A_value,b);
    D optimal;
    D relative_gap;
    D gap = 0;

    L tmp1 = 0;
    L tmp2 = 0;
    L tmp3 = 0;
    std::vector<D> step;
    D eta_sigma;
    D ACDM_tau;
    D ACDM_tol = lp_param.sub_tol;
    Mtranspose(m,n_bar,A_row_ptr,A_col_idx,A_value,AT_row_ptr,AT_col_idx,AT_value);
    Mdiag(n_bar,AT_row_ptr,AT_value,step);
    for (L row = 0; row < n_bar; ++row)
        step[row] +=1;
    std::vector<D> p(n_bar,0);
    std::vector<D> p_sum(n_bar+1,0);

    if(subtype == "exact") {
        time_start = clock();
        AtimesAT(AATplusI,m,n_bar,A_row_ptr,A_col_idx,A_value);
        for (L row = 0; row < m; ++row)
            AATplusI[row * m + row] += 1;
        Matrix_time = (D) (clock() - time_start) / CLOCKS_PER_SEC;
        gsl_matrix_inv(&H.matrix,&Hinv.matrix);
    }else if(subtype == "NU_ACDM"){
        NU_ACDM_param(p,p_sum,ACDM_tau,eta_sigma,n_bar,step);
    }
    D tmp;
    L iter = 0;
    time = 0;
    while((tmp3<n_bar)) {
        ++iter;
        if(subtype == "exact") {
            time_start = clock();
            for (L row = 0; row < n_bar; ++row)
                ve[row] = dual_res_x[row] / rho - y[row];
            for (L row = 0; row < m; ++row)
                we[row] = lambda_E[row] / rho-b[row];
            Axminusb(Aeveminuswe, m, A_row_ptr, A_col_idx, A_value, ve, we);
            Ax(g, m, Inv, Aeveminuswe);
            ATxminusb(x, m, n_bar, A_row_ptr, A_col_idx, A_value, g, ve);
            time += (D) (clock() - time_start) / CLOCKS_PER_SEC;
        }else if(subtype == "NU_ACDM"){
            time_start = clock();
            ATx(ge,m,n_bar, A_row_ptr, A_col_idx, A_value,we);
            for (L row = 0; row < n_bar; ++row)
                ge[row] = - (dual_res[row]/rho - ATb[row] -y[row]);
            time += (D) (clock() - time_start) / CLOCKS_PER_SEC;
            NU_ACDM(x,Atimesx,time,n_bar,m,AT_row_ptr,AT_col_idx,AT_value,step,ge,p,p_sum,ACDM_tol,ACDM_tau,eta_sigma);
        }
        //Step2:
        time_start = clock();
        for (L row = 0; row < nb; ++row) {
            y[row] = x[row] + lambda_x[row]/rho;
            if (y[row] < 0)
                y[row] = 0;
        }
        for (L row = nb; row < n; ++row)
            y[row] = x[row] + lambda_x[row]/rho;

        for (L row = n; row < n_bar; ++row) {
            y[row] = x[row] + lambda_x[row]/rho;
            if (y[row] < 0)
                y[row] = 0;
        }
        if(subtype == "exact")
            Ax(Atimesx,m, A_row_ptr, A_col_idx, A_value, x);
        for (L row = 0; row < m; ++row)
            d[row] = Atimesx[row] - b[row];
        //Axminusb(d, m, A_row_ptr, A_col_idx, A_value, x,b);

        for (L row = 0; row < n_bar; ++row)
            xminusy[row] = x[row] - y[row];

        for (L row = 0; row < m; ++row)
            lambda_E[row] += rho * d[row];
        for (L row = 0; row < n_bar; ++row)
            lambda_x[row] += rho * xminusy[row];
        ATx(dual_res_E,m,n_bar,A_row_ptr,A_col_idx,A_value,lambda_E);
        for (L row = 0; row < n; ++row)
            dual_res_x[row] = lambda_x[row]+c[row];
        for (L row = n; row < n_bar; ++row)
            dual_res_x[row] = lambda_x[row];
        for (int row = 0; row < n_bar; ++row)
            dual_res[row] = dual_res_E[row] + dual_res_x[row];


        tmp1 = 0;
        for (L row = 0; row < m; ++row) {
            if (d[row] * d[row] > tol)
                break;
            ++tmp1;

        }
        tmp2 = 0;
        if (tmp1 == m) {
            for (L row = 0; row < n_bar; ++row) {
                if (xminusy[row] * xminusy[row] > tol)
                    break;
                ++tmp2;
            }
        }
        tmp3 = 0;
        if (tmp2 == n_bar) {
            for (L row = 0; row < n_bar; ++row) {
                if (dual_res[row] * dual_res[row] > tol)
                    break;
                ++tmp3;
            }
        }

       time += (D) (clock() - time_start)/CLOCKS_PER_SEC;
       Print("tmp1",tmp1);
       Print("tmp2",tmp2);
        Print("tmp3",tmp3);
   }

   time += Matrix_time;
   Print("Matrix_time",Matrix_time);
   Print("time",time);
   Print("iter",iter);
    optimal = 0;
       for (L row = 0; row < n; ++row)
           optimal += c[row]*x[row];
       Print("optimal",optimal);
       gap = 0;
       for (L row = 0; row < m; ++row)
           gap += b[row]*lambda_E[row];
       gap += optimal;
       Print("gap",gap);
       relative_gap = sqrt(gap*gap/(optimal*optimal));
       Print("relative_gap",relative_gap);

   return x;
}

template<typename L, typename D>
std::vector<D> LP<L, D>::ADMM2(const D &tol,const std::string &subtype) {
   D rho = lp_param.rho;
   setA2();
   setb2();

   D time = 0;
   D time_start;
   D Matrix_time;
   std::vector<D> z(n, 0);
   std::vector<D> y(n, 0);
   std::vector<D> lambda_I(m_bar, 0);
   std::vector<D> lambda_z(n, 0);
   std::vector<D> AATplusI(m_bar * m_bar, 0);
   std::vector<D> Inv(m_bar * m_bar, 0);
   gsl_matrix_view H = gsl_matrix_view_array(&AATplusI[0], m_bar, m_bar);
   gsl_matrix_view Hinv = gsl_matrix_view_array(&Inv[0], m_bar, m_bar);

   std::vector<D> wi(m_bar, 0);
   std::vector<D> Xi(m_bar, 0);
   std::vector<D> vi(n, 0);
   std::vector<D> d(m_bar, 0);
   std::vector<D> zminusy(n, 0);
   std::vector<D> Atimesz(m_bar,0);
   std::vector<D> Aiviminuswi(m_bar, 0);
   std::vector<D> g(m_bar, 0);
   std::vector<D> gi(n,0);
   std::vector<D> tmp_vec1(m_bar, 0);
   std::vector<D> tmp_vec2(m_bar, 0);
   std::vector<D> dual_res(n,0);
   std::vector<D> dual_res_I(n,0);
   std::vector<D> dual_res_x(n,0);
   D gap = 0;
   D optimal;
   D relative_gap;
   L tmp1 = 0;
   L tmp2 = 0;
   L tmp3 = 0;
   std::vector<D> step;
   D eta_sigma;
   D ACDM_tau;
   D ACDM_tol = lp_param.sub_tol;

   Mtranspose(m_bar,n,A_row_ptr,A_col_idx,A_value,AT_row_ptr,AT_col_idx,AT_value);
   Mdiag(n,AT_row_ptr,AT_value,step);
   for (L row = 0; row < n; ++row)
       step[row] +=1;
   std::vector<D> p(n,0);
   std::vector<D> p_sum(n+1,0);

   if (subtype == "exact") {
       time_start = clock();
       AtimesAT(AATplusI, m_bar, n, A_row_ptr, A_col_idx, A_value);
       for (L row = 0; row < m_bar; ++row)
           AATplusI[row * m_bar + row] += 1;
       Matrix_time = (D) (clock() - time_start) / CLOCKS_PER_SEC;
       gsl_matrix_inv(&H.matrix, &Hinv.matrix);
   }else if(subtype == "NU_ACDM"){
       NU_ACDM_param(p,p_sum,ACDM_tau,eta_sigma,n,step);
   }
   L iter = 0;
   time = 0;
   while((tmp1!=m_bar) || (tmp2!= n)) {
       ++iter;

       //Step1:
       time_start = clock();
       for (L row = 0; row < n; ++row)
           vi[row] = (c[row] + lambda_z[row]) / rho - y[row];
       for (L row = 0; row < m_bar; ++row)
           wi[row] = lambda_I[row] / rho-b[row]+ Xi[row];
       time += (D) (clock() - time_start) / CLOCKS_PER_SEC;
       if(subtype == "exact") {
           time_start = clock();
           Axminusb(Aiviminuswi, m_bar, A_row_ptr, A_col_idx, A_value, vi, wi);
           Ax(g, m_bar, Inv, Aiviminuswi);
           ATxminusb(z, m_bar, n, A_row_ptr, A_col_idx, A_value, g, vi);
           time += (D) (clock() - time_start) / CLOCKS_PER_SEC;
       }else if(subtype == "NU_ACDM"){
           time_start = clock();
           ATx(gi,m_bar,n, A_row_ptr, A_col_idx, A_value,wi);
           for (L row = 0; row < n; ++row)
               gi[row] = - gi[row] - vi[row];
           time += (D) (clock() - time_start) / CLOCKS_PER_SEC;
           NU_ACDM(z,Atimesz,time,n,m_bar,AT_row_ptr,AT_col_idx,AT_value,step,gi,p,p_sum,ACDM_tol,ACDM_tau,eta_sigma);

       }
       //Step2:
       time_start = clock();
       for (L row = 0; row < nb; ++row) {
           y[row] = z[row] + lambda_z[row]/rho;
           if (y[row] < 0)
               y[row] = 0;
       }
       for (L row = nb; row < n; ++row)
           y[row] = z[row] + lambda_z[row]/rho;
       if(subtype == "exact")
           Ax(Atimesz,m_bar, A_row_ptr, A_col_idx, A_value, z);
       for (L row = 0; row < m; ++row)
           d[row] = Atimesz[row] - b[row];
       //Axminusb(d, m_bar, A_row_ptr, A_col_idx, A_value, z,b);
       for (L row = 0; row < m_bar; ++row) {
           tmp_vec1[row] = d[row] + lambda_I[row] / rho;
           if(tmp_vec1[row] < 0)
               Xi[row] = -tmp_vec1[row];
           else
               Xi[row] = 0;
           tmp_vec2[row] = d[row] + Xi[row];
       }

       for (L row = 0; row < n; ++row)
           zminusy[row] = z[row] - y[row];

       for (L row = 0; row < m_bar; ++row)
           lambda_I[row] += rho * tmp_vec2[row];

       for (L row = 0; row < n; ++row)
           lambda_z[row] += rho * zminusy[row];

       tmp1 = 0;
       for (L row = 0; row < m_bar; ++row) {
           if (tmp_vec2[row] * tmp_vec2[row] > tol)
               break;
           ++tmp1;
       }
       tmp2 = 0;
       if (tmp1 == m_bar) {
           for (L row = 0; row < n; ++row) {
               if (zminusy[row] * zminusy[row] > tol)
                   break;
               ++tmp2;
           }
       }
        time += (D) (clock() - time_start)/CLOCKS_PER_SEC;
        Print("tmp1",tmp1);
        Print("tmp2",tmp2);
        Print("tmp3",tmp3);

    }
    time += Matrix_time;
    Print("Matrix_time",Matrix_time);
    Print("time",time);
    Print("iter",iter);
    optimal = 0;
    for (L row = 0; row < n; ++row)
        optimal += c[row]*z[row];
    gap = 0;
    for (L row = 0; row < n; ++row)
        gap += b[row]*
    Print("optimal");
    return z;
}

/*
template<typename L, typename D>
std::vector<D> LP<L, D>::SSNPPA(const D &tol) {
    std::vector<D> x(n,0);
    setAi();
    setbi();

    std::vector<D> wi(mi,0);
    std::vector<D> lambda_I(mi,0);
    std::vector<D> we(me,0);
    std::vector<D> lambda_E(me,0);
    std::vector<D> J;
    std::vector<D> AJ_row_ptr;
    std::vector<D> AJ_col_idx;
    std::vector<D> AJ_value;
    std::vector<D> wJ;
    std::vector<D> Grad;
    D eta = 1;
    L sizeJ;

    Axminusb(we,me, Ae_row_ptr,Ae_col_idx,Ae_value,x,be);
    for (L row = 0; row < me; ++row)
        we[row] += lambda_E[row];

    for (L row = 0; row < mi; ++row)
        wi[row] += lambda_I[row];
    Axminusb(wi,mi, Ai_row_ptr,Ai_col_idx,Ai_value,x,bi);
    for (L row = 0; row < mi; ++row)
        wi[row] += lambda_I[row];
    J.clear();
    for (L row = 0; row < me; ++row)
        J.push_back(row);
    wJ = we;
    for (L row = 0; row < mi; ++row) {
        if(wi[row] > 0) {
            J.push_back(row + me);
            wJ.push_back(wi[row]);
        }
    }
    sizeJ = J.size();
    subMat(AJ_row_ptr,AJ_col_idx,AJ_value,A_col_idx,A_row_ptr,A_value,sizeJ,J);




    return x;
}
 */




/*
template<typename L, typename D>
std::vector<D> LP<L, D>::ADMM(const D &tol) {
    //parameter for ADMM
    D rho = lp_param.rho;
    D tau = lp_param.tau;
    D sigma = tau*rho;

    //parameter for inner iterations
    D sub_tol = lp_param.sub_tol;


    //primal variables
    std::vector<D> x(n,0);
    std::vector<D> y(n,0);
    std::vector<D> z(n,0);
    std::vector<D> Xi(mi,0);
    std::vector<D> xminusy(n,0);
    std::vector<D> zminusy(n,0);

    //dual variables
    std::vector<D> lambda_I(mi,0);
    std::vector<D> lambda_E(me,0);
    std::vector<D> lambda_x(n,0);
    std::vector<D> lambda_z(n,0);

    //variables for equality
    std::vector<D> we(me,0);
    std::vector<D> ve(n,0);
    std::vector<D> de(me,0);
    std::vector<D> tmpe1(me,0);
    std::vector<D> tmpe2(me,0);
    std::vector<D> tmpe3(n,0);
    std::vector<D> tmpi1(mi,0);
    std::vector<D> tmpi2(mi,0);
    std::vector<D> He(me*me,0);
    std::vector<D> Hi(mi*mi,0);
    std::vector<D> InvHe(me*me,0);
    std::vector<D> InvHi(mi*mi,0);
    gsl_matrix_view gsl_He = gsl_matrix_view_array(&He[0],me,me);
    gsl_matrix_view gsl_InvHe = gsl_matrix_view_array(&InvHe[0],me,me);
    gsl_matrix_view gsl_Hi = gsl_matrix_view_array(&Hi[0],mi,mi);
    gsl_matrix_view gsl_InvHi = gsl_matrix_view_array(&InvHi[0],mi,mi);

    AtimesAT(He,me,n,Ae_row_ptr,Ae_col_idx,Ae_value);
    for (L row = 0; row < me; ++row)
        He[row*me+row] +=1;
    gsl_matrix_inv(&gsl_He.matrix, &gsl_InvHe.matrix);


    AtimesAT(Hi,mi,n,Ai_row_ptr,Ai_col_idx,Ai_value);
    for (L row = 0; row < mi; ++row)
        Hi[row*mi+row] +=1;
    gsl_matrix_inv(&gsl_Hi.matrix, &gsl_InvHi.matrix);
    //Print(mi,InvHi);

    //variables for inequality
    std::vector<D> wi(mi,0);
    std::vector<D> vi(n,0);
    std::vector<D> di(mi,0);

    L tmp1=0,tmp2=0,tmp3=0,tmp4=0;
    D gap;
    while(tmp4 < n){
        //step 1: update x and z

        //update x exactly: x =-ve + AeT*(InvHe*(Ae*ve-we));
        for (L row = 0; row < me; ++row)
            we[row] = lambda_E[row]/rho - be[row];
        for (L row = 0; row < n; ++row)
            ve[row] = (lambda_x[row] + c[row])/rho - y[row];

        Axminusb(tmpe1,me,Ae_row_ptr,Ae_col_idx,Ae_value,ve,we);
        Ax(tmpe2,me,InvHe,tmpe1);
        ATxminusb(x,me,n,Ae_row_ptr,Ae_col_idx,Ae_value,tmpe2,ve);

        for (L row = 0; row < mi; ++row)
            wi[row] = di[row] + lambda_I[row]/rho + Xi[row];
        for (L row = 0; row < n; ++row)
            vi[row] = zminusy[row] + (lambda_z[row] +c[row])/rho;

        Axminusb(tmpe1,mi,Ai_row_ptr,Ai_col_idx,Ai_value,vi,wi);
        Ax(tmpe2,mi,InvHi,tmpe1);
        ATxminusb(x,mi,n,Ai_row_ptr,Ai_col_idx,Ai_value,tmpe2,vi);


        //step 2: update y = 0.5( x + z - (lambda_x + lambda_z)/rho )
        for (L row = 0; row < n; ++row)
            y[row] = 0.5*(x[row] + z[row] +(lambda_x[row] + lambda_z[row])/rho);
        for (L row = 0; row < nb; ++row) {
            if(y[row] < 0)
                y[row] = 0;
        }

        //step 3: update lambda
        //de = Ae*x - be; di = Ai*z -bi;
        Axminusb(de,me,Ae_row_ptr,Ae_col_idx,Ae_value,x,be);
        Axminusb(di,mi,Ai_row_ptr,Ai_col_idx,Ai_value,z,bi);
        for (L row = 0; row < mi; ++row) {
            tmpi1[row] = di[row] + lambda_I[row] / rho;
            if(tmpi1[row] < 0)
                Xi[row] = -tmpi1[row];
            else
                Xi[row] = 0;
            tmpi2[row] = di[row] + Xi[row];
        }
        for (L row = 0; row < n; ++row) {
            xminusy[row] = x[row] - y[row];
            zminusy[row] = z[row] - y[row];
        }
        for (L row = 0; row < me; ++row)
            lambda_E[row] += sigma*de[row];
        for (L row = 0; row < mi; ++row)
            lambda_I[row] = lambda_I[row] + sigma*tmpi2[row];
        for (L row = 0; row < n; ++row) {
            lambda_x[row] += sigma*xminusy[row];
            lambda_z[row] += sigma*zminusy[row];
        }

        tmp1 = 0;
        tmp2 = 0;
        tmp3 = 0;
        tmp4 = 0;
        gap = 0;
        for (L row = 0; row < n; ++row)
            gap += c[row] * (z[row] + x[row]);
        for (int row = 0; row < me; ++row)
            gap += be[row]*lambda_E[row];
        for (L row = 0; row < mi; ++row)
            gap += bi[row] * lambda_I[row];
        gap *= 0.5;
        Print("gap",gap);
        if(gap*gap >tol) {
            while (de[tmp1] * de[tmp1] <= tol)
                ++tmp1;
            if (tmp1 == me) {
                while (xminusy[tmp2] * xminusy[tmp2] <= tol)
                    ++tmp2;
            }
            if (tmp2 == n) {
                while ((tmpi2[tmp3] * tmpi2[tmp3] <= tol))
                    ++tmp3;
            }
            if (tmp3 == mi) {
                while (zminusy[tmp4] * zminusy[tmp4] <= tol)
                    ++tmp4;
            }
        }
    }

    return x;
}
 */


#endif //SDNA_GENINFO_H
