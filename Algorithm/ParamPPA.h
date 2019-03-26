//
// Created by Lu Meng on 2019/3/2.
//

#ifndef ALGORITHM_PARAMPPA_H
#define ALGORITHM_PARAMPPA_H

#include "../Problem/ProblemHeader.h"
template<typename L, typename D>
inline void exp_update(D& varepsilon, const D& varepsilon_initial, const D& exp_idx, const L& iter);

template<typename L, typename D>
inline void inv_update(D& varepsilon, const D& varepsilon_initial, const D& inv_idx, const L& iter);

template<typename D>
class ParamPPA {
public:
    //exp_idx < 1; cur_idx >1;
    //initial varepsilon and the parameters for updating varepsilon
    D varepsilon;
    D exp_idx_varepsilon = 0.9;
    D inv_idx_varepsilon = 1.1;

    //initial delta and the parameters for updating delta
    D delta;
    D exp_idx_delta = 0.9;
    D inv_idx_delta = 1.1;

    /* theta: initial guess for theta of max(Hoffman constant, growth condition parameter)
     * sigma = alpha * theta: the parameter for PPA, i.e, P = (I+sigma*T)^{-1}
     * exp_idx_sigma: we also can say exp_idx_theta, every time we double theta, then we also double sigma
     */
    D theta = 1.;
    D alpha;
    //sigma = alpha * theta;
    D sigma;
    D exp_idx_theta = 2;

    //Notica that beta_multiple in (0,1)
    D beta;
    D beta_multiple;
    ParamPPA(const D& alpha = 1., const D&beta_multiple = 0.5){
        this->alpha = alpha;
        this->beta_multiple = beta_multiple;
        /* If we set alpha = 1., beta_multiple = 0.5
         * then beta = 0.1733 and delta = 0.0471
         */
        D alpha_square = alpha * alpha;
        beta = beta_multiple * log(sqrt(1+alpha_square));
        delta = (exp(-beta)- 1./sqrt(1+alpha_square))/(2+exp(-beta));
        varepsilon = delta;
    }
    ~ParamPPA(){}

    void printparam(){
        Print("varepsilon_initial", varepsilon);
        Print("exp_idx_varepsilon", exp_idx_varepsilon);
        Print("inv_idx_varepsilon", inv_idx_varepsilon);

        Print("delta_initial", delta);
        Print("exp_idx_delta", exp_idx_delta);
        Print("inv_idx_delta", inv_idx_delta);

        Print("theta_intial", theta);
        Print("exp_idx_theta", exp_idx_theta);
        Print("alpha", alpha);
        Print("sigma_initial", sigma);
        Print("beta_multiple", beta_multiple);
        Print("beta", beta);
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


#endif //ALGORITHM_PARAMPPA_H
