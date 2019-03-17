//
// Created by Lu Meng on 2018/10/24.
//

#ifndef USE_QUADPROGFUN_H
#define USE_QUADPROGFUN_H

#include "useheader.h"
#include "Print.h"
#include "../QuadProg/src/QuadProg++.hh"

template <typename L, typename D>
quadprogpp::Matrix<D> VtoM_qp(const L &n, const L &d, const std::vector<D> &v){
    quadprogpp::Matrix<D> G;
    G.resize(n,d);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < d; ++j) {
            G[i][j] = v[i*d+j];
        }
    }
    return G;
};

template <typename L, typename D>
quadprogpp::Vector<D> VtoV_qp(const L &d, const std::vector<D> &v){
    quadprogpp::Vector<D> V;
    V.resize(d);
    for (int i = 0; i < d; ++i) {
        V[i] = v[i];
    }
    return V;
};


template <typename L, typename D>
void Print(const L &n, const L &d, const quadprogpp::Matrix<D> &G){
    for (L i = 0; i < n; ++i) {
        for (L j = 0; j < d; ++j)
            std::cout << G[i][j] << " ";
        Print();
    }
};

template <typename L, typename D>
void Print(const L &d, const quadprogpp::Vector<D> &G){
    for (L i = 0; i < d; ++i) {
        std::cout << G[i] << " ";
    }
    Print();
};



#endif //USE_QUADPROGFUN_H
