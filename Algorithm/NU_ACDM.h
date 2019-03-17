#ifndef ALGORITHM_NU_ACDM_H
#define ALGORITHM_NU_ACDM_H

#include "AbsAlgorithm.h"
#include "AbsAlgorithm.h"

template<typename L, typename D>
class NU_ACDM:public AbsAlgorithm<L, D>{
private:
    D beta = 0;
    D gamma, S_gamma, tau, eta;
    L m,n,K;
    std::vector<D> y;
    std::vector<D> z;
    std::vector<D> w_y;
    std::vector<D> w_z;
    std::vector<D> cord_grad_sum;
    std::vector<D> delta_y;
    std::vector<D> delta_z;
    std::vector<L> cord;
    std::vector<D> w_before;
    std::vector<D> x_before;
    std::vector<D> p;
    std::vector<D> p_sum;
    D fun_value_before;
    D precision;
    bool stop;
    gsl_rng *gsl_rng_r;
public:
    NU_ACDM(ProbData<L, D>* const data_inst, const L& loss_type);
    ~NU_ACDM() {
        gsl_rng_free(gsl_rng_r);
    }
    void sub_initial();
    void sub_step();
    void step();
    void solver(const std::vector<D>& x_initial, const L& max_iter, const D& precision,const L& blocksize);
    void left_update();
    void result_print();


};
template<typename L, typename D>
NU_ACDM<L, D>::NU_ACDM(ProbData<L, D> *const data_inst, const L &loss_type):
        AbsAlgorithm<L, D>(data_inst, loss_type) {
    gsl_rng_env_setup();
    gsl_rng_r = gsl_rng_alloc(gsl_rng_default);
    m = data_inst->m;
    n = data_inst->n;
    y = std::vector<D>(n);
    z = std::vector<D>(n);
    w_y = std::vector<D>(m);
    w_z = std::vector<D>(m);
    w_before = std::vector<D>(m);
    x_before = std::vector<D>(n);
    AbsAlgorithm<L,D>::grad = std::vector<D>(n);
    AbsAlgorithm<L,D>::w = std::vector<D>(m);
    AbsAlgorithm<L,D>::x = std::vector<D>(n);

    gamma = 0.5*(1-beta);
    p = std::vector<D>(n);
    for (L row = 0; row < n; ++row)
        p[row] = pow(AbsAlgorithm<L,D>::Lip[row], gamma);
    S_gamma = 0;
    for (L row = 0; row < n; ++row)
        S_gamma += p[row];
    for (L row = 0; row < n; ++row)
        p[row] /= S_gamma;
    for (L row = 0; row < n; ++row)
        p_sum[row + 1] = p_sum[row] + p[row];

    tau = 2. / (1 + sqrt(4 * S_gamma * S_gamma/ alpha + 1));
    eta = scale / (tau * S_gamma*S_gamma);

}

template<typename L, typename D>
void NU_ACDM<L, D>::sub_initial() {
    theta = tau/(D)n;
    for (L row = 0; row < n; ++row){
        y[row] = 0.;
        z[row] =  AbsAlgorithm<L,D>::x[row];
    }
    for (L row = 0; row < m; ++row){
        w_y[row] = 0;
        w_z[row] =  AbsAlgorithm<L,D>::w[row];
    }
}

template<typename L, typename D>
void NU_ACDM<L, D>::sub_step() {
    D tmp = n* theta/(D)tau;
    D theta_square = theta*theta;
    AbsAlgorithm<L,D>::loss->cord_grad_sum_update(cord_grad_sum,AbsAlgorithm<L,D>::data,w_y,y,w_z,z,theta_square,tau,cord);
    //Print("cord_grad_sum",cord_grad_sum);
    for (L row = 0; row < tau ; ++row) {
        delta_z[row] = -cord_grad_sum[row];
        delta_z[row] /= (tmp*AbsAlgorithm<L,D>::v[cord[row]]);
    }
    //Print("t",delta_z);
    D tmp_u = -(1- tmp)/theta_square;
    //Print("z",z);
    for (L row = 0; row < tau ; ++row) {
        if((cord[row] < AbsAlgorithm<L,D>::data->nplus) && (z[cord[row]] + delta_z[row] < 0)){
            z[cord[row]] = 0;
            delta_z[row] = - z[cord[row]];
        }else{
            z[cord[row]] += delta_z[row];
        }
        delta_y[row] = tmp_u * delta_z[row];
        y[cord[row]] += delta_y[row];

    }
    //Print("cord",cord);
    //Print("t",delta_z);
    AbsAlgorithm<L,D>::loss->w_update(w_z,AbsAlgorithm<L,D>::data,delta_z,tau,cord);
    AbsAlgorithm<L,D>::loss->w_update(w_y,AbsAlgorithm<L,D>::data,delta_y,tau,cord);
    theta = 0.5*(sqrt(theta_square*theta_square + 4*theta_square) - theta_square);
    /*
    Print("cord",cord);
    Print("delta_z",delta_z);
    Print("y", y);
    Print("z", z);
    Print("w_y", w_y);
    Print("w_z", w_z);
     */

}

template<typename L, typename D>
void NU_ACDM<L, D>::step() {
    sub_initial();
    for (L i = 0; i < K; ++i) {
        SampleGen(cord,tau,n,gsl_rng_r);
        sub_step();
    }
    D tmp = theta*theta/(1-theta);
    for (L row = 0; row < n; ++row)
        AbsAlgorithm<L,D>::x[row] = tmp*y[row] + z[row];
    for (L row = 0; row < m; ++row)
        AbsAlgorithm<L,D>::w[row] = tmp*w_y[row] + w_z[row];
    AbsAlgorithm<L,D>::loss->fun_value_compute(AbsAlgorithm<L,D>::fun_value,AbsAlgorithm<L,D>::data,
                                               AbsAlgorithm<L,D>::w, AbsAlgorithm<L,D>::x);
    if( AbsAlgorithm<L,D>::fun_value > fun_value_before){
        for (L row = 0; row < n; ++row)
            AbsAlgorithm<L,D>::x[row] = x_before[row];
        for (L row = 0; row < m; ++row)
            AbsAlgorithm<L,D>::w[row] = w_before[row];
    }else{
        for (L row = 0; row < n; ++row)
            x_before[row] =  AbsAlgorithm<L,D>::x[row];
        for (L row = 0; row < m; ++row)
            w_before[row] =  AbsAlgorithm<L,D>::w[row];
        fun_value_before =  AbsAlgorithm<L,D>::fun_value;
        AbsAlgorithm<L,D>::loss->grad_update(AbsAlgorithm<L,D>::grad,AbsAlgorithm<L,D>::data, AbsAlgorithm<L,D>::w,
                                             AbsAlgorithm<L,D>::x);
        left_update();
    }
    //Print("grad_norm:",AbsAlgorithm<L, D>::grad_norm);
    ++AbsAlgorithm<L,D>::iter;
}


template<typename L, typename D>
void NU_ACDM<L, D>::solver(const std::vector<D>& x_initial,const L& max_iter, const D& precision_value,const L& blocksize){
    tau = blocksize;
    cord_grad_sum = std::vector<D>(tau);
    delta_y = std::vector<D>(tau);
    delta_z = std::vector<D>(tau);
    cord = std::vector<L>(tau);
    AbsAlgorithm<L, D>::set_v_mu(tau);
    //Print("v:",AbsAlgorithm<L,D>::v);
    K = ceil(2 * exp(1) * AbsAlgorithm<L,D>::data->n * (sqrt(1+ 1./AbsAlgorithm<L,D>::mu) -1)/(D)tau + 1);
    //Print("K:",K);
    precision = precision_value;
    D tmp;
    for (L row = 0; row < n; ++row) {
        tmp = x_initial[row];
        AbsAlgorithm<L,D>::x[row] = tmp;
        x_before[row] = tmp;
    }
    AbsAlgorithm<L,D>::loss->w_initial( AbsAlgorithm<L,D>::w,AbsAlgorithm<L,D>::data, AbsAlgorithm<L,D>::x,
                                        AbsAlgorithm<L,D>::loss->b);
    for (L row = 0; row < m; ++row)
        w_before[row] =  AbsAlgorithm<L,D>::w[row];

    AbsAlgorithm<L,D>::loss->fun_value_compute(fun_value_before,AbsAlgorithm<L,D>::data, AbsAlgorithm<L,D>::w,
                                               AbsAlgorithm<L,D>::x);
    AbsAlgorithm<L,D>::loss->grad_update(AbsAlgorithm<L,D>::grad,AbsAlgorithm<L,D>::data, AbsAlgorithm<L,D>::w,
                                         AbsAlgorithm<L,D>::x);
    left_update();
    //Print("w", AbsAlgorithm<L,D>::w);
    //Print("grad",AbsAlgorithm<L,D>::grad);
    AbsAlgorithm<L,D>::iter = 0;
    stop = false;
    while((AbsAlgorithm<L, D>::grad_norm > precision) && ( AbsAlgorithm<L,D>::iter < max_iter))
        step();
}

template<typename L, typename D>
void NU_ACDM<L, D>::left_update() {
    D tmp;
    D left = 0;
    for (L row = 0; row < n; ++row) {
        tmp = AbsAlgorithm<L, D>::grad[row];
        left +=  tmp*tmp;
    }
    AbsAlgorithm<L, D>::grad_norm = sqrt(left);
}

template<typename L, typename D>
void NU_ACDM<L, D>::result_print() {
    Print("NU_ACDM_outer_iter:", AbsAlgorithm<L, D>::iter);
    Print("NU_ACDM_grad_norm:",AbsAlgorithm<L, D>::grad_norm);
    Print("NU_ACDM_fun_value:",AbsAlgorithm<L, D>::fun_value);
    Print("NU_ACDM_x:",AbsAlgorithm<L, D>::x);
    Print("NU_ACDM_w:",AbsAlgorithm<L, D>::w);
    Print("NU_ACDM_grad:",AbsAlgorithm<L, D>::grad);
}



#endif //ALGORITHM_NU_ACDM_H
