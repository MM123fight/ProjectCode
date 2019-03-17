//
// Created by Lu Meng on 2019/2/20.
//

#ifndef ALGORITHM_ADAPPPA_H
#define ALGORITHM_ADAPPPA_H

#include "AbsAlgorithmSub.h"
#include "AbsAlgorithm.h"
#include "ParamAdapPPA.h"

template<typename L, typename D>
class AdapPPA:public AbsAlgorithmSub<L, D>{
private:
    L m,mi,n,tau;
    L sub_iter,K;
    L sub_max_iter = 1000;
    //primal_dual__gap,residual;
    D gap, res, res_before;
    D sigma,theta,C,CK_param;
    D sub_precision,sub_alpha,right;
    D varepsilon_stage1, varepsilon_stage2, varepsilon_stage3, epsilon;

    bool stop;
    std::vector<D> sub_b;
    std::vector<D> sub_d;
    std::vector<D> x_before;
    std::vector<D> x_diff_over_sigma;
    std::vector<D> lambda_over_sigma;
    std::vector<D> lambda_over_sigma_before;
    ParamAdapPPA<D>* param = NULL;
public:
    AdapPPA(ProbData<L, D>* const data_inst, const unsigned int& sub_method_type,
            const unsigned int& sub_loss_type_value);

    virtual ~AdapPPA(){
        delete param;
    }
    void sub_initial();
    void sub_step();
    void step();
    void solver(const std::vector<D>& x_initial, const std::vector<D> &lambda_initial,
                const L& max_iter, const D& precision, const L& blocksize);
    void right_update();
    void res_update();
    void res_initial();
};

template<typename L, typename D>
AdapPPA<L, D>::AdapPPA(ProbData<L, D>* const data_inst, const unsigned int& sub_method_type,
                       const unsigned int& sub_loss_type_value):
        AbsAlgorithmSub<L,D>(data_inst,sub_method_type, sub_loss_type_value){
    param = new ParamAdapPPA<D>();
    m = data_inst->m;
    mi = data_inst->mi;
    n = data_inst->n;
    sub_b = std::vector<D>(m);
    sub_d = std::vector<D>(n);
    x_before = std::vector<D>(n);
    x_diff_over_sigma = std::vector<D>(n);
    lambda_over_sigma = std::vector<D>(m);
    lambda_over_sigma_before = std::vector<D>(m);

    AbsAlgorithmSub<L, D>::x = std::vector<D>(n);
    AbsAlgorithmSub<L, D>::w = std::vector<D>(m);
    AbsAlgorithmSub<L, D>::lambda = std::vector<D>(m);
    AbsAlgorithmSub<L, D>::res_dual = std::vector<D>(n);
    D tmp1 = 0;
    D tmp2 = 0;
    D tmp;
    for (L row = 0; row < n; ++row) {
        tmp = data_inst->c[row];
        tmp1 += tmp*tmp;
    }
    for (L row = 0; row < m; ++row) {
        tmp =  data_inst->b[row];
        tmp2 += tmp*tmp;
    }
    tmp = 2. * std::max<D>(tmp1,tmp2);
    std::vector<D> v(n);
    AbsAlgorithmSub<L,D>::set_loss();
    AbsAlgorithmSub<L, D>::loss->v_intial_set(v,data_inst,1);
    for (L row = 0; row < n; ++row)
        tmp += v[row];
    CK_param = sqrt(tmp);
    //Print("v",v);
    //data_inst->dataprint();
    //Print("CK_param",CK_param);
}

//fix sub_b, sub_d, epsilon_k, delta_k,
template<typename L, typename D>
void AdapPPA<L, D>::sub_initial() {
    sigma = param->alpha * theta;
    sub_alpha = 1./(sigma*sigma);
    C = theta * CK_param * res;
    K = ceil(log(theta * CK_param)/param->beta);
    inv_update(varepsilon_stage1, param->varepsilon_initial,param->inv_idx,AbsAlgorithmSub<L, D>::iter);
    //Print("C",C);
    Print("K",K);
    Print("sigma",sigma);
    Print("varepsilon_stage1",varepsilon_stage1);
    // exp_update(varepsilon_stage1, param->varepsilon_initial,param->exp_idx,AbsAlgorithmSub<L, D>::iter);
}

//PPA_step
template<typename L, typename D>
void AdapPPA<L, D>::sub_step() {
    Print("sub_iter", sub_iter);
    inv_update<L,D>(varepsilon_stage3, varepsilon_stage2, param->inv_idx,sub_iter);
    Print("varepsilon_stage3",varepsilon_stage3);
    sub_precision = sub_alpha * varepsilon_stage3;
    //Print("sub_precision",sub_precision);
    for (L row = 0; row < m; ++row)
        sub_b[row] = AbsAlgorithmSub<L,D>::data->b[row] - lambda_over_sigma[row];
    for (L row = 0; row < n; ++row)
        sub_d[row] = AbsAlgorithmSub<L,D>::x[row] -sigma*AbsAlgorithmSub<L,D>::data->c[row];
    //Print("sub_b",sub_b);
    //Print("sub_d",sub_d);
    AbsAlgorithmSub<L,D>::sub_method->set_loss(sub_alpha,sub_b,sub_d);
    AbsAlgorithmSub<L,D>::sub_method->solver(AbsAlgorithm<L,D>::x, sub_max_iter, sub_precision, tau);
    right_update();
    while(AbsAlgorithmSub<L,D>::sub_method->grad_norm > right) {
        AbsAlgorithmSub<L,D>::sub_method->step();
         if(AbsAlgorithmSub<L,D>::sub_method->grad_norm <= sub_precision)
             right_update();
    }
    //update x
    for (L row = 0; row < n; ++row) {
        x_diff_over_sigma[row] = (AbsAlgorithmSub<L,D>::sub_method->x[row] - AbsAlgorithmSub<L,D>::x[row])/sigma;
        AbsAlgorithmSub<L, D>::x[row] = AbsAlgorithmSub<L, D>::sub_method->x[row];
    }
    //update lambda_over_sigma
    D tmp;
    for (L row = 0; row < mi; ++row){
        tmp = AbsAlgorithmSub<L,D>::sub_method->w[row];
        AbsAlgorithmSub<L,D>::w[row] = tmp - lambda_over_sigma[row];
        if(tmp < 0)
            tmp = 0;
        lambda_over_sigma[row] = tmp;
    }
    for (L row = mi; row < m; ++row) {
        tmp = AbsAlgorithmSub<L,D>::sub_method->w[row];
        AbsAlgorithmSub<L,D>::w[row] = tmp - lambda_over_sigma[row];
        lambda_over_sigma[row] = tmp;
    }
    ++sub_iter;
    //Print("lambda_over_sigma",lambda_over_sigma);
}

template<typename L, typename D>
void AdapPPA<L, D>::step() {
    sub_initial();
    L t = 0;
    D res_bound = C;
    D res_before = res;
    //Print("residual_before", res);
    while( (res <= res_bound)&&(res<=res_before)&&(res>epsilon) ) {
        Print("t",t);
        res_before = res;
        for (L row = 0; row < n; ++row)
            x_before[row] = AbsAlgorithmSub<L,D>::x[row];
        for (L row = 0; row < m; ++row)
            lambda_over_sigma_before[row] = lambda_over_sigma[row];
        inv_update(varepsilon_stage2, varepsilon_stage1,param->inv_idx,t);
        Print("varepsilon_stage2",varepsilon_stage2);
        sub_iter = 0;
        while (sub_iter < K) {
            sub_step();
        }
        res_update();
        ++t;
        res_bound = C*exp(-param->beta * K * t);
        Print("residual", res);
        //Print("x",AbsAlgorithmSub<L,D>::x);
        //Print("lambda",AbsAlgorithmSub<L,D>::lambda);
    }

    if(res <= epsilon) {
        Print("exit with residual no more than epsilon");
        exit(0);
    }
    else if(res > res_before) {
        res = res_before;
        for (L row = 0; row < n; ++row)
            AbsAlgorithmSub<L,D>::x[row] = x_before[row];
        for (L row = 0; row < m; ++row)
            lambda_over_sigma[row] = lambda_over_sigma_before[row];
    }
    theta *= param->multiple;
    ++AbsAlgorithmSub<L,D>::iter;
}

template<typename L, typename D>
void AdapPPA<L, D>::solver(const std::vector<D> &x_initial, const std::vector<D> &lambda_initial,
                           const L &max_iter, const D &precision, const L &blocksize) {
    //param->printparam();
    //AbsAlgorithmSub<L,D>::data->dataprint();
    tau = blocksize;
    epsilon = precision;
    theta = param->theta_initial;
    for (L row = 0; row < n; ++row)
        AbsAlgorithmSub<L,D>::x[row] = x_initial[row];
    for (L row = 0; row < m; ++row) {
        lambda_over_sigma[row] = lambda_initial[row]/(param->alpha * theta);
    }
    AbsAlgorithmSub<L,D>::loss->w_initial(AbsAlgorithmSub<L,D>::w,AbsAlgorithmSub<L,D>::data, x_initial,
                                          AbsAlgorithmSub<L,D>::data->b);
    AbsAlgorithmSub<L,D>::loss->res_dual_initial(AbsAlgorithmSub<L,D>::res_dual,AbsAlgorithmSub<L,D>::data,
                                                 lambda_initial,AbsAlgorithmSub<L,D>::data->c);
    res_initial();
    Print("residual_intial", res);
    AbsAlgorithmSub<L,D>::iter = 0;
    while((res > epsilon)&&(AbsAlgorithmSub<L,D>::iter < max_iter)) {
        step();
    }
    AbsAlgorithmSub<L,D>::primal_dual_gap = gap;
}

template<typename L, typename D>
void AdapPPA<L, D>::right_update() {
    right = 0;
    D tmp;
    for (L row = 0; row < n; ++row) {
        tmp = (AbsAlgorithmSub<L,D>::sub_method->x[row] - AbsAlgorithmSub<L,D>::x[row])/sigma;
        right += tmp*tmp;
    }
    //AbsAlgorithmSub<L,D>::sub_method->result_print();

    for (L row = 0; row < mi; ++row) {
        tmp = AbsAlgorithmSub<L,D>::sub_method->w[row];
        if(tmp < 0) {
            tmp = 0;
        }
        tmp -= lambda_over_sigma[row];
        right += tmp * tmp;
    }
    for (L row = mi; row < m; ++row) {
        tmp = AbsAlgorithmSub<L,D>::sub_method->w[row] - lambda_over_sigma[row];
        right += tmp * tmp;
    }
    right = sqrt(right);
    right *= (param->delta/sigma);
    Print("sub_outer_iter",AbsAlgorithmSub<L,D>::sub_method->iter);
    Print("sub_grad_norm",AbsAlgorithmSub<L,D>::sub_method->grad_norm);
    Print("right", right);
}




template<typename L, typename D>
void AdapPPA<L, D>::res_update() {

    //primal_dual_gap
    gap = 0;
    for (L row = 0; row < m; ++row)
        gap += AbsAlgorithmSub<L,D>::data->b[row]*lambda_over_sigma[row];
    gap *= sigma;
    for (L row = 0; row < n; ++row)
        gap += AbsAlgorithmSub<L,D>::data->c[row] * AbsAlgorithmSub<L,D>::x[row];
    res = gap*gap;
    D tmp;
    //primal_res
    for (L row = 0; row < mi; ++row){
        tmp = AbsAlgorithmSub<L,D>::w[row];
        if(tmp > 0)
            res += tmp*tmp;
    }
    for (L row = mi; row < m; ++row) {
        tmp = AbsAlgorithmSub<L,D>::w[row];
        res += tmp*tmp;
    }
    //dual_res
    for (L row = 0; row < n; ++row) {
        tmp = sigma * AbsAlgorithmSub<L, D>::sub_method->grad[row] - x_diff_over_sigma[row];
        res += tmp * tmp;
    }
    res = sqrt(res);
    AbsAlgorithmSub<L, D>::residual = res;
}

template<typename L, typename D>
void AdapPPA<L, D>::res_initial() {
    AbsAlgorithmSub<L, D>::sub_method->grad.clear();
    AbsAlgorithmSub<L, D>::sub_method->grad.resize(n,0);
    x_diff_over_sigma.clear();
    x_diff_over_sigma.resize(n,0);

    res_update();
    D tmp;
    res = res*res;
    for (L row = 0; row < n; ++row) {
        tmp = AbsAlgorithmSub<L, D>::res_dual[row];
        res += tmp * tmp;
    }
    res = sqrt(res);
    AbsAlgorithmSub<L, D>::residual = res;
}


#endif //ALGORITHM_ADAPPPA_H
