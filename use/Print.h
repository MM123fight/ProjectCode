//
// Created by Lu Meng on 2018/10/23.
//

#ifndef USE_PRINT_H
#define USE_PRINT_H

#include "stdafx.h"

int Prec = 22;

//Print nothing but a empty line
void  Print();
//Print values or string
template <typename D>
void  Print(const D &v);
//Print vector
template <typename D>
void  Print(const std::vector<D> &v);
//Print two segments of values or string
void Print(const std::string &u, const std::string &v);
//Print string:vector
template <typename D>
void Print(const std::string &u, const std::vector<D> &v);
//Print matrix
template <typename L, typename D>
void Print(const L &n, const L &d, const std::vector<L> &row_ptr, const std::vector<L> &col_idx,
           const std::vector<D> &value, const std::string mode="Full");


void  Print() {
    std::cout << std::endl;
}

template <typename D>
void  Print(const D &v){
    std::cout.precision(Prec);
    std::cout << v << std::endl;
}


void Print(const std::string &u, const std::string &v){
    std::cout << u << " " << v << std::endl;
}

template <typename D>
void Print(const std::string &u, const std::vector<D> &v){
    std::cout << u << " ";
    Print(v);
}

template <typename D>
void  Print(const std::vector<D> &v){
    std::cout.precision(Prec);
    for(int i = 0; i < v.size(); ++i)
        std::cout << v[i] << " ";
    std::cout << std::endl;
}


template <typename L, typename D>
void Print(const L &n, const L &d, const std::vector<L> &row_ptr, const std::vector<L> &col_idx,
           const std::vector<D> &value, const std::string mode){
    if(mode == "Full"){
        std::vector<D> tmp = std::vector<D>(d);
        for (L row = 0; row < n; ++row) {
            tmp.clear();
            tmp.resize(d,0);
            for (L col = row_ptr[row]; col < row_ptr[row+1]; ++col)
                tmp[col_idx[col]] = value[col];
            Print(tmp);
        }
    }else if(mode == "Sparse"){
        for (L row = 0; row < n; ++row) {
            for (L col = row_ptr[row]; col < row_ptr[row+1]; ++col)
                std::cout << col_idx[col]+1 << ":" << value[col] << " ";
            Print();
        }
    }
    else{
        Print("The mode is not available");
    }
}

template<typename D>
void Print(const int &tau, const std::vector<D> &value) {
    for (int row = 0; row < tau; ++row) {
        for (int col = 0; col < tau; ++col) {
            std::cout << value[row*tau+col] << " ";
        }
        Print();
    }
}

template<typename D>
void Print(const std::string name ,const  D z) {
    std::cout << name << ":" << z << std::endl;
}




#endif
