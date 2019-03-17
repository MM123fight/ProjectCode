#ifndef ALGORITHM_APPROX_H
#define ALGORITHM_APPROX_H

#include "AbsAlgorithm.h"
#include "PCDM.h"

template<typename L, typename D>
class APPROX:public AbsAlgorithm<L, D>{
protected:
    L m,n,tau,K;
    std::vector<D> y;
    std::vector<D> z;
    std::vector<D> w_y;
    std::vector<D> w_z;
    std::vector<D> cord_grad_sum;
    std::vector<D> delta_y;
    std::vector<D> delta_z;
    std::vector<L> cord;

    D theta, theta_square;
    D precision;
    gsl_rng *gsl_rng_r;
    std::string method_type;
public:
    APPROX(ProbData<L, D>* const data_inst, const unsigned int& loss_type_value, const unsigned int& reg_type_value,
           const std::string& data_file_path="");
    virtual ~APPROX() {
        gsl_rng_free(gsl_rng_r);
    }
    void sub_step();
    void step();
    void sub_initial();
    void restart_step();
    void solver(const std::vector<D>& x_initial, const L& blocksize, const L& max_iter,
                const D& precision_value, const unsigned int& stop_type = 1, const bool& result_iflog_value = false,
                const std::string& method_type = "");
    void initial_print();
    void result_print();


};
template<typename L, typename D>
APPROX<L, D>::APPROX(ProbData<L, D> *const data_inst, const unsigned int &loss_type_value,const unsigned int &reg_type_value,
                     const std::string& data_file_path_value):
        AbsAlgorithm<L, D>(data_inst, loss_type_value,reg_type_value,data_file_path_value) {
    gsl_rng_env_setup();
    gsl_rng_r = gsl_rng_alloc(gsl_rng_default);
    m = data_inst->m;
    n = data_inst->n;
    y = std::vector<D>(n);
    z = std::vector<D>(n);
    w_y = std::vector<D>(m);
    w_z = std::vector<D>(m);
    AbsAlgorithm<L,D>::grad = std::vector<D>(n);
    AbsAlgorithm<L,D>::w = std::vector<D>(m);
    AbsAlgorithm<L,D>::x = std::vector<D>(n);
}


template<typename L, typename D>
void APPROX<L, D>::sub_step() {
    D tmp = n* theta/(D)tau;
    AbsAlgorithm<L,D>::loss->cord_grad_sum_update(cord_grad_sum,AbsAlgorithm<L,D>::data,w_y,y,w_z,z,theta_square,tau,cord);
    D tmp_u = -(1- tmp)/theta_square;
    AbsAlgorithm<L, D>::reg_cord_cal(z,delta_z,tmp,cord_grad_sum,tau,cord);
    for (L row = 0; row < tau ; ++row){
        delta_y[row] = tmp_u * delta_z[row];
        y[cord[row]] += delta_y[row];
    }
    AbsAlgorithm<L,D>::loss->w_update(w_z,AbsAlgorithm<L,D>::data,delta_z,tau,cord);
    AbsAlgorithm<L,D>::loss->w_update(w_y,AbsAlgorithm<L,D>::data,delta_y,tau,cord);
    theta = 0.5*(sqrt(theta_square*theta_square + 4*theta_square) - theta_square);
    theta_square = theta * theta;
}

template<typename L, typename D>
void APPROX<L, D>::step() {
    D start = clock();
    for (L i = 0; i < K; ++i) {
        SampleGen(cord,tau,n,gsl_rng_r);
        sub_step();
    }
    D tmp = theta_square/(1-theta_square);
    std::vector<D> sum(n);
    std::vector<D> sum_w(m);
    std::vector<D> sum_grad(n);
    D sum_grad_norm;
    for (L row = 0; row < n; ++row) {
        AbsAlgorithm<L, D>::x[row] = tmp * y[row] + z[row];
        sum[row] = theta_square * y[row] + z[row];
    }
    for (L row = 0; row < m; ++row) {
        AbsAlgorithm<L, D>::w[row] = tmp * w_y[row] + w_z[row];
        sum_w[row] = theta_square * w_y[row] + w_z[row];
    }
    AbsAlgorithm<L, D>::grad_compute(sum_grad,sum_w,sum);
    sum_grad_norm = 0;
    L max_cord = 0;
    D max_cord_grad = sum_grad[0];
    for (int j = 0; j < n; ++j) {
        if(sum_grad[j] < max_cord_grad) {
            max_cord_grad = sum_grad[j];
            max_cord = j;
        }
        sum_grad_norm += sum_grad[j] * sum_grad[j];
    }
    AbsAlgorithm<L, D>::left_value_compute();
    ++AbsAlgorithm<L,D>::iter;
    AbsAlgorithm<L, D>::time += (D)(clock() - start)/CLOCKS_PER_SEC;
    if(AbsAlgorithm<L, D>::result_iflog == true)
        AbsAlgorithm<L, D>::result_log();
}


template<typename L, typename D>
void APPROX<L, D>::sub_initial() {
    theta = tau/(D)n;
    theta_square = theta * theta;
    y.clear();
    y.resize(n,0);
    z.assign(AbsAlgorithm<L,D>::x.begin(),AbsAlgorithm<L,D>::x.end());
    w_y.clear();
    w_y.resize(m,0);
    w_z.assign(AbsAlgorithm<L,D>::w.begin(),AbsAlgorithm<L, D>::w.end());
}

template<typename L, typename D>
void APPROX<L, D>::restart_step(){
    sub_initial();
    step();
}

template<typename L, typename D>
void APPROX<L, D>::solver(const std::vector<D>& x_initial, const L& blocksize, const L& max_iter,
                          const D& precision_value, const unsigned int& stop_type_value, const bool& result_iflog_value,
                          const std::string& method_type_value){
    method_type = method_type_value;
    Print("method_type",method_type);
    AbsAlgorithm<L, D>::check_blocksize(blocksize);
    tau = blocksize;
    AbsAlgorithm<L, D>::result_iflog = result_iflog_value;
    if(AbsAlgorithm<L, D>::result_iflog == true) {
        AbsAlgorithm<L, D>::full_result_path = AbsAlgorithm<L, D>::data->data_file_path + "/result/APPROX_bs" + to_string(tau)
                                               + "_alpha" + to_string(AbsAlgorithm<L, D>::loss->alpha);
        AbsAlgorithm<L, D>::logfile_open();
    }

    AbsAlgorithm<L, D>::stop_type = stop_type_value;
    precision = precision_value;
    AbsAlgorithm<L, D>::set_v_mu(tau);
    AbsAlgorithm<L, D>::Krcd_compute(K,tau);


    cord_grad_sum = std::vector<D>(tau);
    delta_y = std::vector<D>(tau);
    delta_z = std::vector<D>(tau);
    cord = std::vector<L>(tau);
    AbsAlgorithm<L,D>::x.assign(x_initial.begin(), x_initial.end());
    AbsAlgorithm<L,D>::primal_initial_set();

    //Print("left",AbsAlgorithm<L, D>::left);
    AbsAlgorithm<L, D>::iter = 0;
    AbsAlgorithm<L, D>::time = 0;
    if(AbsAlgorithm<L, D>::result_iflog == true)
        AbsAlgorithm<L, D>::result_log();
    if(method_type == "restart"){
        while((AbsAlgorithm<L, D>::left > precision) && ( AbsAlgorithm<L,D>::iter < max_iter)) {
            restart_step();
        }
    }else{
        sub_initial();
        while((AbsAlgorithm<L, D>::left > precision) && ( AbsAlgorithm<L,D>::iter < max_iter)) {
            step();
        }
    }
    AbsAlgorithm<L, D>::logFile.close();
}


template<typename L, typename D>
void APPROX<L, D>::initial_print() {
    if(method_type == "restart"){
        Print("ResAPPROX_initial_fun_value:",AbsAlgorithm<L, D>::initial_fun_value);
        Print("ResAPPROX_initial_grad_norm:",AbsAlgorithm<L, D>::initial_grad_norm);
    } else {
        Print("APPROX_initial_fun_value:", AbsAlgorithm<L, D>::initial_fun_value);
        Print("APPROX_initial_grad_norm:", AbsAlgorithm<L, D>::initial_grad_norm);
    }
}

template<typename L, typename D>
void APPROX<L, D>::result_print() {
    AbsAlgorithm<L, D>::fun_value_compute(AbsAlgorithm<L, D>::fun_value, AbsAlgorithm<L, D>::w, AbsAlgorithm<L, D>::x);
    AbsAlgorithm<L, D>::grad_compute(AbsAlgorithm<L, D>::grad, AbsAlgorithm<L, D>::w, AbsAlgorithm<L, D>::x);
    AbsAlgorithm<L, D>::grad_norm_update(AbsAlgorithm<L, D>::grad_norm);
    if(method_type == "restart"){
        Print("ResAPPROX_outer_iter:", AbsAlgorithm<L, D>::iter);
        Print("ResAPPROX_time:", AbsAlgorithm<L, D>::time);
        Print("ResAPPROX_fun_value:", AbsAlgorithm<L, D>::fun_value);
        Print("ResAPPROX_grad_norm:", AbsAlgorithm<L, D>::grad_norm);
        Print("ResAPPROX_x:", AbsAlgorithm<L, D>::x);
        Print("ResAPPROX_w:", AbsAlgorithm<L, D>::w);
        Print("ResAPPROX_grad:", AbsAlgorithm<L, D>::grad);

    }else {
        Print("APPROX_outer_iter:", AbsAlgorithm<L, D>::iter);
        Print("APPROX_time:", AbsAlgorithm<L, D>::time);
        Print("APPROX_fun_value:", AbsAlgorithm<L, D>::fun_value);
        Print("APPROX_grad_norm:", AbsAlgorithm<L, D>::grad_norm);
        Print("APPROX_x:", AbsAlgorithm<L, D>::x);
        Print("APPROX_w:", AbsAlgorithm<L, D>::w);
        Print("APPROX_grad:", AbsAlgorithm<L, D>::grad);
    }
}



#endif //ALGORITHM_RESAPPROX_H
