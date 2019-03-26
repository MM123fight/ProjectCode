//
// Created by Lu Meng on 2019/3/19.
//

#include <iostream>
#include <vector>
#include <math.h>
#include "APPROX.h"

#include "../use/useheader.h"

template <typename L>
void Sample(std::vector<L>& cord, L tau, const L& n);
template <typename L, typename D>
inline void norm_square(D& sqaure, std::vector<D>& vec){
    sqaure = 0.;
    D tmp;
    for (L row = 0; row < vec.size(); ++row) {
        tmp = vec[row];
        sqaure += tmp*tmp;
    }
}


int main()
{
    int tau= 2;
    int n = 100;
    std::vector<int> A(4,1);
    A[0] = 1;
    A[1] = 2;
    A[2] = 3;
    A[3] = 4;
    std::vector<int> C = A;
    std::vector<int> S(4);

    double time;
    int square;
    int sample_times = 6000000;
    double start;


    time = 0.;

    time = 0.;
    double theta = 1.;
    time =0.;
    start = clock();

    int mi = tau/2;
    int tmp;

    time = 0.;
    start = clock();
    for (int times = 0; times < sample_times*10; ++times) {
        MtimesM_s2(C,A);
    }
    time += (double)(clock()-start)/CLOCKS_PER_SEC;
    std::cout << "time:" << time << std::endl;
    Print(C);

    time = 0;
    C= A;
    start = clock();
    for (int times = 0; times < sample_times*10; ++times) {
        MtimesM_s2(S,C,A);
        C = S;
    }
    time += (double)(clock()-start)/CLOCKS_PER_SEC;
    std::cout << "time:" << time << std::endl;

    Print(C);








   // Print("time",time);

    return 0;
}

template <typename L>
void Sample(std::vector<L>& cord, L tau, const L& n){
    L j = 0;
    for (L i = 0; i < n; i++) {
        if (rand() % (n - i) < tau) {
            cord[j] = i;
            tau--;
            j++;
        }
    }
}


