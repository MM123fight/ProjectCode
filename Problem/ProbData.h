#ifndef PROBDATA_H
#define PROBDATA_H

#include "../use/useheader.h"
#include "LPuse.h"

//First construct LP to LS problem
/* LP problem: A is a n*m matrix         can be transformed into a Least square problem with
 * A:n*d, b:R^n, x:R^d, c:R^d                    (size nb) (size m)  (size nf)
 *                                       Ai_LS =   Ai[nb]      0       Ai[nf]      x_LS = x_[nb]       bi_LS = bi
 * min <c,x>                             Ae_LS =   Ae[nb]      0       Ae[nf]        lambda_[mi]       be_LS = be
 *                                                    0       A^T        0           lambda_[me]               -c
 * Ae*x = be                                       c[nb]^T     b       c[nf]^T            x_[nf]                0
 * Ai*x <= bi                            mplus = nb + mi, nplus = nb+mi, m_LS = m + n + 1; n_LS = n+m
 * x_[nb] >= 0
 *
 *
 *                                       AT_LS =  [nb]  AT_nb    [n]   0     c
 *                                                [m]     0      [m]   A     b
 *                                                [nf]   AT_nf
 */


template <typename L, typename D>
class ProbData{
public:

    L mi,me,m;
    L nb,nf,n;
    L nplus,mplus;

    //Ai: mi*n
    std::vector<L> Ai_row_ptr = std::vector<L>();
    std::vector<L> Ai_col_idx = std::vector<L>();
    std::vector<D> Ai_value = std::vector<D>();


    // Ae: me*n
    std::vector<L> Ae_row_ptr = std::vector<L>();
    std::vector<L> Ae_col_idx = std::vector<L>();
    std::vector<D> Ae_value = std::vector<D>();


    // A = [Ai;Ae], m*n
    std::vector<L> A_row_ptr = std::vector<L>();
    std::vector<L> A_col_idx = std::vector<L>();
    std::vector<D> A_value = std::vector<D>();

    // AT: n*m
    std::vector<L> AT_row_ptr = std::vector<L>();
    std::vector<L> AT_col_idx = std::vector<L>();
    std::vector<D> AT_value = std::vector<D>();

    // bi: R^mi, be:R^me
    // b = [bi;be], R^m
    std::vector<D> bi = std::vector<D>();
    std::vector<D> be = std::vector<D>();
    std::vector<D> b = std::vector<D>();
    //c: R^n
    std::vector<D> c = std::vector<D>();
    std::string data_type;
    std::string data_file_path;

    ProbData(const std::string &root_path_name, const std::string &data_file_name, std::string data_type = "LP");
    virtual ~ProbData(){}
    void data_print();

    void set_mplus(const D& mplus) {
        if (mplus <= m) {
            this->mplus = mplus;
        } else {
            Print("The mplus exceed m! Please reset mplus");
            exit(0);
        }
    }
    void set_nplus(const D& nplus){
        if (nplus <= n) {
            this->nplus = nplus;
        } else {
            Print("The nplus exceed n! Please reset nplus");
            exit(0);
        }
    }
};

template <typename L, typename D>
ProbData<L, D>::ProbData(const std::string &root_path_name, const std::string &data_file_name, std::string data_type) {
    this->data_type = data_type;
    data_file_path = root_path_name + "/"+ data_file_name;
    Print("Loading Data...");
    if((data_type == "LS")||(data_type == "LS_equ")) {
        MatrixRead(data_file_path, m, n, A_row_ptr, A_col_idx, A_value, b);
        nplus = n;
        mplus = m;
        if(data_type == "LS")
            Mtranspose(AT_row_ptr, AT_col_idx, AT_value, m, n, A_row_ptr, A_col_idx, A_value);
    }else if((data_type == "LP")||(data_type == "LPtoLP_inequ")||(data_type == "LPtoLS")||(data_type == "LPtoLS_equ")){
        //Read Ai and Ae
        readMeta(data_file_path+"/meta", mi, me, nb, nf);
        n = nb+nf;
        readMat(data_file_path+"/A", mi, n,Ai_row_ptr, Ai_col_idx, Ai_value);
        readMat(data_file_path+"/Aeq", me, n,Ae_row_ptr, Ae_col_idx, Ae_value);
        VectorRead(data_file_path+"/b",mi,bi);
        VectorRead(data_file_path+"/beq",me,be);
        VectorRead(data_file_path+"/c",n,c);
        //A=[Ai;Ae];
        if(data_type == "LPtoLP_inequ"){
            std::vector<D> zeros(nb,0);
            std::vector<D> ones(nb,-1);
            std::vector<L> ones_idx(nb);
            for (L row = 0; row < nb; ++row)
                ones_idx[row] = row;
            Ai_col_idx.insert(Ai_col_idx.begin(),ones_idx.begin(),ones_idx.end());
            Ai_value.insert(Ai_value.begin(),ones.begin(),ones.end());
            bi.insert(bi.begin(),zeros.begin(),zeros.end());
            for (L row = 0; row < mi+1; ++row)
                Ai_row_ptr[row] += nb;
            Ai_row_ptr.insert(Ai_row_ptr.begin(),ones_idx.begin(),ones_idx.end());
            mi += nb;
            nb = 0;
            nf = n;
            /*
            Print("Ai_row_ptr",Ai_row_ptr);
            Print("Ai_col_idx",Ai_col_idx);
            Print("Ai_value",Ai_value);
             */
        }
        m = mi+me;
        nplus = nb;
        mplus = mi;

        A_row_ptr = Ai_row_ptr;
        A_col_idx = Ai_col_idx;
        A_value = Ai_value;
        b = bi;
        A_col_idx.insert(A_col_idx.end(), Ae_col_idx.begin(),Ae_col_idx.end());
        A_value.insert(A_value.end(), Ae_value.begin(), Ae_value.end());
        L tmp = A_row_ptr.back();
        for (L row = 1; row < me+1; ++row)
            A_row_ptr.push_back(Ae_row_ptr[row]+tmp);
        b.insert(b.end(),be.begin(),be.end());
        Mtranspose(AT_row_ptr,AT_col_idx,AT_value,m,n,A_row_ptr,A_col_idx,A_value);

        if((data_type == "LPtoLS")||(data_type == "LPtoLS_equ")){


            /*can be transformed into a Least square problem with
                         [nb]  [m]   [nf]
               Ai_LS = Ai[nb]   0    Ai[nf]      x_LS = x_[nb]           bi_LS = bi
               Ae_LS =    0    -AT     0               lambda_[mi]            c
                       Ae[nb]   0    Ae[nf]            lambda_[me]       be_LS = be
                       c[nb]    b     c[nf]             x_[nf]                0
             */
            for (L row = 0; row < m; ++row) {
                for (L col = A_row_ptr[row]; col < A_row_ptr[row+1]; ++col) {
                    if(A_col_idx[col] >= nb)
                        A_col_idx[col] += m;
                }
            }

            for (L row = 0; row < n; ++row) {
                for (L col = AT_row_ptr[row]; col < AT_row_ptr[row+1]; ++col) {
                   AT_col_idx[col] += nb;
                   AT_value[col] = -AT_value[col];
                }
            }

            L tmp = A_row_ptr[mi];
            Print("tmp",tmp);
            for(L row = 1; row < n+1; ++row)
                AT_row_ptr[row] += tmp;
            A_row_ptr.insert(A_row_ptr.begin() + mi+1,AT_row_ptr.begin()+1,AT_row_ptr.end());
            A_col_idx.insert(A_col_idx.begin() + tmp,AT_col_idx.begin(),AT_col_idx.end());
            A_value.insert(A_value.begin() + tmp,AT_value.begin(),AT_value.end());
            L tmp2 = A_row_ptr[mi+n];
            for (L row = mi+n+1; row < m + n+ 1; ++row) {
                A_row_ptr[row] += tmp2 - tmp;
            }

            std::vector<L> full_idx(m+n);
            for (L row = 0; row < m + n; ++row)
                full_idx[row] = row;
            A_value.insert(A_value.end(),c.begin(),c.begin()+nb);
            A_value.insert(A_value.end(),b.begin(),b.end());
            A_value.insert(A_value.end(),c.begin()+nb,c.end());
            A_col_idx.insert(A_col_idx.end(),full_idx.begin(),full_idx.end());
            A_row_ptr.push_back(A_row_ptr.back()+m+n);
            b.insert(b.begin() + mi,c.begin(),c.end());
            b.push_back(0);
            n = m + n;
            m = n + 1;
            mplus = nb + mi;
            nplus = nb + mi;
            //Print("A_row_ptr", A_row_ptr);
            //Print("A_col_idx", A_col_idx);
            //Print("A_value", A_value);
            if(data_type == "LPtoLS") {
                Mtranspose(AT_row_ptr, AT_col_idx, AT_value, m, n, A_row_ptr, A_col_idx, A_value);
                Print("check");
                std::string LPtoSVM_path = data_file_path+"/LPtoSVM";
                SvmWrite(LPtoSVM_path,A_row_ptr,A_col_idx,A_value,b);
                Print("check");
            }
        }
    }else{
        Print("There is no such data type");
        exit(0);
    }
    if((data_type == "LS_equ")||(data_type == "LPtoLS_equ")){
        for (L row = 0; row < A_col_idx.size(); ++row)
            A_col_idx[row] += mplus;
        L idx;
        for (L row = 0; row < mplus; ++row) {
            idx = A_row_ptr[row]+row;
            A_value.insert(A_value.begin()+idx,1.);
            A_col_idx.insert(A_col_idx.begin()+idx,row);
        }
        for (L row = 1; row <= mplus; ++row)
            A_row_ptr[row] += row;
        for (L row = mplus+1; row < m+1; ++row)
            A_row_ptr[row] += mplus;
        n = mplus + n;
        nplus += mplus;
        mplus = 0;
        Mtranspose(AT_row_ptr, AT_col_idx, AT_value, m, n, A_row_ptr, A_col_idx, A_value);
    }
    Print("Finish loading data from:", data_file_path);

}

template <typename L, typename D>
void ProbData<L, D>::data_print(){
    if(data_type == "LP"||"LPtoLP_inequ") {
        std::cout << "m=" << m << " mi=" << mi << " me=" << me << " mplus=" << mplus << std::endl;
        std::cout << "n=" << n << " nb=" << nb << " nf=" << nf << " nplus=" << nplus << std::endl;
        Print("c",c);
    }else{
        std::cout << "m=" << m << " n=" << n << " mplus=" << mplus <<  " nplus=" << nplus << std::endl;
    }
    Print("b",b);
    Print("Matrix A:");
    Print(m,n,A_row_ptr,A_col_idx,A_value);
    Print("Matrix AT:");
    Print(n,m,AT_row_ptr,AT_col_idx,AT_value);
    if(data_type == "LP"||"LPtoLP_inequ") {
        Print("Matrix Ai:");
        Print(mi,n,Ai_row_ptr,Ai_col_idx,Ai_value);
        Print("Matrix Ae:");
        Print(me,n,Ae_row_ptr,Ae_col_idx,Ae_value);
    }
}

#endif //PROBDATA_H
