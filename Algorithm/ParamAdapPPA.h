//
// Created by Lu Meng on 2019/3/2.
//

#ifndef ALGORITHM_PARAMADAPPPA_H
#define ALGORITHM_PARAMADAPPPA_H

#include "../Problem/ProblemHeader.h"
template<typename L, typename D>
inline void exp_update(D& varepsilon, const D& varepsilon_initial, const D& exp_idx, const L& iter);

template<typename L, typename D>
inline void inv_update(D& varepsilon, const D& varepsilon_initial, const D& inv_idx, const L& iter);

template<typename D>
class ParamAdapPPA {
public:
    //exp_idx < 1; cur_idx >1;
    D varepsilon_initial = 1.;
    D exp_idx = 0.9;
    D inv_idx = 1.1;

    D theta_initial = 1.;
    D multiple = 2.;

    D alpha = 1.;
    D beta;
    D delta;

    ParamAdapPPA(){
        D alpha_square = alpha * alpha;
        beta = 0.5*log(sqrt(1+alpha_square));
        delta = (exp(-beta)- 1./sqrt(1+alpha_square))/(2+exp(-beta));
    }
    ~ParamAdapPPA(){}

    void printparam(){
        Print("exp_idx", exp_idx);
        Print("inv_idx", inv_idx);
        Print("theta_intial", theta_initial);
        Print("varepsilon_initial", varepsilon_initial);
        Print("alpha", alpha);
        Print("beta", beta);
        Print("delta", delta);
    }
};

template<typename L, typename D>
inline void exp_update(D& varepsilon, const D& varepsilon_initial, const D& exp_idx, const L& iter){
    varepsilon = varepsilon_initial * pow(exp_idx,iter);
}

template<typename L, typename D>
inline void inv_update(D& varepsilon, const D& varepsilon_initial, const D& inv_idx, const L& iter){
    varepsilon = varepsilon_initial/pow(iter+1,inv_idx);
}


#endif //ALGORITHM_PARAMADAPPPA_H
