//
// Created by Lu Meng on 2018/10/23.
//

#ifndef USE_SAMPLE_H
#define USE_SAMPLE_H

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_linalg.h>
#include "stdafx.h"
#include "Print.h"

//p_sum = [0,p1,p1+p2,...,p1+...+p_(n-1),1]
template<typename L, typename D>
inline L RandSample(const L& n, const std::vector<D> &p_sum,gsl_rng *gsl_rng_r){
    //generate [0,1)
    D tmp = gsl_rng_uniform(gsl_rng_r);
    L check_point_start = 0;
    L check_point_end = n;
    L check_point = (check_point_start+check_point_end)/2+1;
    while(check_point != check_point_start) {
        if(tmp == p_sum[check_point]) {
            break;
        } else if(tmp > p_sum[check_point]) {
            check_point_start = check_point;
        } else{
            check_point_end = check_point;
        }
        check_point = (check_point_start+check_point_end)/2;
    }
    return check_point+1;
};

template<typename L>
inline std::vector<L> SampleGen(const L &max, const int &tau, gsl_rng *gsl_rng_r){
    std::vector<L> S = std::vector<L>(tau);
    // Calculate a tau-nice sampling
    if (tau < max) {
        for (int i = 0; i < tau; ++i) {
            bool done = true;
            do {
                done = true;
                S[i] = gsl_rng_uniform_int(gsl_rng_r, max);
                for (int j = 0; j < i; ++j) {
                    if (S[i] == S[j]) {
                        done = false;
                        break;
                    }
                }

            } while (!done);
        }
    } else {
        for (int i = 0; i < tau; i++)
            S[i] = i;
    }
    return S;
}





template<typename L>
inline void SampleGen(std::vector<L> &cord, const L& blocksize, const L& max, gsl_rng *gsl_rng_r) {
    if (blocksize == 1){
        cord[0] = gsl_rng_uniform_int(gsl_rng_r, max);
    }else if(blocksize < max) {
        for (int i = 0; i < blocksize; ++i) {
            bool done = true;
            do {
                done = true;
                cord[i] = gsl_rng_uniform_int(gsl_rng_r, max);
                for (int j = 0; j < i; ++j) {
                    if (cord[i] == cord[j]) {
                        done = false;
                        break;
                    }
                }

            } while (!done);
        }
    } else {
        for (int i = 0; i < blocksize; i++)
            cord[i] = i;
    }
}

#endif //USE_SAMPLE_H
