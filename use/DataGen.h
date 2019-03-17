//
// Created by Lu Meng on 2018/10/23.
//
//Use matlab to create your data with the package libsvm;
#ifndef USE_DATAGEN_H
#define USE_DATAGEN_H

#include "Sample.h"
#include "usefun.h"
#include "Print.h"
template <typename L, typename D>
void DataGen(L m, L n, D nnz_percent, std::string data_dir)
{
    std::vector< std::vector<L> > nnz_index(m); //Index of matrix A;
    std::vector< std::vector<D> > nnz_value(m); //Value of matrix A;
    std::vector<L> row_nnz_num(m); //Number of each row of A
    std::vector<D> d(m,0); // vector b;

    /*for(L i = 0; i < m; i++)
        row_nnz_num[i] = (rand()%maxnum_nnz)+1;*/
    for(L i = 0; i < m; i++)
        row_nnz_num[i] = nnz_percent * n;

    //Assume x_opt[i] = 1, for all i = 0, ...
    L real_n = 0;
    gsl_rng_env_setup();
    gsl_rng *gsl_rng_r = gsl_rng_alloc(gsl_rng_default);
    for(L i = 0; i < m; i++)
    {
        nnz_index[i]=SampleGen(n, row_nnz_num[i],gsl_rng_r);
        nnz_value[i].resize(row_nnz_num[i]);
        for(L j = 0; j < row_nnz_num[i]; j++)
        {
            //Generate value randomly from [-1,1];
            nnz_value[i][j] = 2*(rand()%1000)/1000.0-1;
            d[i] += nnz_value[i][j];
        }
        //Calculate the real n;
        if(nnz_index[i][row_nnz_num[i]-1] > real_n)
            real_n = nnz_index[i][row_nnz_num[i]-1];
    }
    ++real_n;
    gsl_rng_free(gsl_rng_r);

    //make b with ||b|| =1;
    D d_sum = 0 ;
    for(L i = 0; i < m; i++)
        d_sum = d_sum + d[i]*d[i];
    d_sum = sqrt(d_sum);
    for(L i = 0; i < m; i++)
    {
        d[i] /= d_sum;
        for(L j = 0; j < row_nnz_num[i]; j++)
            nnz_value[i][j] /= d_sum;
    }

    // Output m_s*n_s matrix
    std::string m_s = to_string(m);
    std::string n_s = to_string(real_n);

    std::string pathnameLS = data_dir + "/data/LS_"+ m_s + "_" + n_s +".txt";
    std::ofstream fout;
    fout.open(pathnameLS.c_str());
    for(L i = 0; i < m ; i++)
    {
        fout << d[i] << " ";
        for(L j = 0; j < row_nnz_num[i]; j++)
            fout << nnz_index[i][j]+1 << ":" << nnz_value[i][j] << " ";
        fout << std::endl;
    }
    fout.close();
    Print("hello");

}

template <typename L, typename D>
void DataGenVector(L n, std::string data_dir){
    std::string n_s = to_string(n);
    std::string pathnameLP_c= data_dir + "/data/c_" + n_s +".txt";
    std::vector<D> c(n);
    D norm = 0;
    std::ofstream fout;
    fout.open(pathnameLP_c.c_str());
    for(L i = 0; i < n ; ++i) {
        c[i] = 2 * (rand() % 1000) / 1000.0 - 1;
        norm += c[i]*c[i];
    }
    norm = sqrt(norm);
    for (L i = 0; i < n; ++i) {
        c[i] /= norm;
        fout << c[i] << std::endl;
    }
    fout.close();
};


#endif //USE_DATAGEN_H
