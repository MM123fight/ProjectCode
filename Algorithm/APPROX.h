#ifndef ALGORITHM_APPROX_H
#define ALGORITHM_APPROX_H

#include "AbsAlgorithm.h"
#include "PCDM.h"


template<typename L, typename D>
class APPROX:public AbsAlgorithm<L, D>{
protected:
    L m,n,tau;
    std::vector<D> y;
    std::vector<D> z;
    std::vector<D> w_y;
    std::vector<D> w_z;
    std::vector<D> cord_grad_sum;
    std::vector<D> delta_y;
    std::vector<D> delta_z;
    std::vector<L> cord;
    std::vector<D> zeros_m;
    std::vector<D> zeros_n;

    D theta, theta_square;
    gsl_rng *gsl_rng_r;
public:
    APPROX(ProbData<L, D>* const data_inst, const unsigned int& loss_type, const unsigned int& reg_type,
           const std::string& data_file_path="");
    virtual ~APPROX() {
        gsl_rng_free(gsl_rng_r);
    }
    void sub_step();
    void step();
    void sub_initial();
    void restart_step();
    void solver(const std::vector<D>& x_initial, const L& blocksize, const L& max_iter,
                const D& precision, const unsigned int& stop_type = 1, const bool& result_iflog = false,
                const std::string& method_type = "");
    void result_print();

};
template<typename L, typename D>
APPROX<L, D>::APPROX(ProbData<L, D> *const data_inst, const unsigned int &loss_type,const unsigned int &reg_type,
                     const std::string& data_file_path):
        AbsAlgorithm<L, D>(data_inst, loss_type,reg_type,data_file_path) {
    m = data_inst->m;
    n = data_inst->n;
    zeros_m = std::vector<D>(m,0);
    zeros_n = std::vector<D>(n,0);
    y = std::vector<D>(n);
    w_y = std::vector<D>(m);
    z = std::vector<D>(n);
    w_z = std::vector<D>(m);
    AbsAlgorithm<L,D>::grad = std::vector<D>(n);
    AbsAlgorithm<L,D>::w = std::vector<D>(m);
    AbsAlgorithm<L,D>::x = std::vector<D>(n);
}


template<typename L, typename D>
void APPROX<L, D>::sub_step() {
    //y: u in paper; z: \tilde{z} in paper
    D tmp = n* theta/(D)tau;
    D tmp_u = -(1- tmp)/theta_square;

    //calculate coordinate grad of f at point theta^2*y +z
    AbsAlgorithm<L,D>::loss->cord_grad_sum_update(cord_grad_sum,AbsAlgorithm<L,D>::data,w_y,y,w_z,z,theta_square,tau,cord);
    //calculate z and delta_z: t in paper;
    AbsAlgorithm<L, D>::reg_cord_cal(z,delta_z,tmp,cord_grad_sum,tau,cord);
    //update y
    for (L row = 0; row < tau ; ++row){
        delta_y[row] = tmp_u * delta_z[row];
        y[cord[row]] += delta_y[row];
    }
    //update w_z and w_y
    AbsAlgorithm<L,D>::loss->w_update(w_z,AbsAlgorithm<L,D>::data,delta_z,tau,cord);
    AbsAlgorithm<L,D>::loss->w_update(w_y,AbsAlgorithm<L,D>::data,delta_y,tau,cord);
    //update theta
    theta = 0.5*(sqrt(theta_square*theta_square + 4*theta_square) - theta_square);
    theta_square = theta * theta;
}

template<typename L, typename D>
void APPROX<L, D>::step() {
    for (L i = 0; i < AbsAlgorithm<L, D>::K; ++i) {
        //Generate Sample
        SampleGen(cord, tau, n, gsl_rng_r);
        sub_step();
    }
    //calculate norm of grad at point x: need full dimension
    D tmp = theta_square/(1-theta_square);
    for (L row = 0; row < n; ++row) {
        AbsAlgorithm<L, D>::x[row] = tmp * y[row] + z[row];
    }
    for (L row = 0; row < m; ++row) {
        AbsAlgorithm<L, D>::w[row] = tmp * w_y[row] + w_z[row];
    }
    AbsAlgorithm<L, D>::left_value_compute();
    ++AbsAlgorithm<L,D>::iter;
}


template<typename L, typename D>
void APPROX<L, D>::sub_initial(){
    theta = tau/(D)n;
    theta_square = theta * theta;
    y = zeros_n;
    w_y = zeros_m;
    z = AbsAlgorithm<L,D>::x;
    w_z = AbsAlgorithm<L,D>::w;
}

template<typename L, typename D>
void APPROX<L, D>::restart_step(){
    sub_initial();
    step();
}

template<typename L, typename D>
void APPROX<L, D>::solver(const std::vector<D>& x_initial, const L& blocksize, const L& max_iter,
                          const D& precision, const unsigned int& stop_type, const bool& result_iflog,
                          const std::string& method_type){
    gsl_rng_env_setup();
    gsl_rng_r = gsl_rng_alloc(gsl_rng_default);
    AbsAlgorithm<L, D>::method_type = method_type;
    //Print("APPROX_method_type",AbsAlgorithm<L, D>::method_type);
    AbsAlgorithm<L, D>::check_blocksize(blocksize);
    AbsAlgorithm<L, D>::blocksize = blocksize;
    tau = blocksize;
    AbsAlgorithm<L, D>::result_iflog = result_iflog;
    if(AbsAlgorithm<L, D>::result_iflog == true) {
        if(method_type == "restart") {
            AbsAlgorithm<L, D>::full_result_path = AbsAlgorithm<L, D>::data->data_file_path + "/result/ResAPPROX_bs" +
                                                   to_string(tau) + "_alpha" + to_string(AbsAlgorithm<L, D>::loss->alpha);
        }else{
            AbsAlgorithm<L, D>::full_result_path = AbsAlgorithm<L, D>::data->data_file_path + "/result/APPROX_bs" +
                                                   to_string(tau) + "_alpha" + to_string(AbsAlgorithm<L, D>::loss->alpha);
        }
        AbsAlgorithm<L, D>::logfile_open();
    }
    AbsAlgorithm<L, D>::stop_type = stop_type;
    AbsAlgorithm<L, D>::precision = precision;
    AbsAlgorithm<L, D>::max_iter = max_iter;
    AbsAlgorithm<L, D>::set_v_mu(tau);
    AbsAlgorithm<L, D>::Krcd_compute(AbsAlgorithm<L, D>::K,tau);
    //Print("K", AbsAlgorithm<L, D>::K);

    cord_grad_sum = std::vector<D>(tau);
    delta_y = std::vector<D>(tau);
    delta_z = std::vector<D>(tau);
    cord = std::vector<L>(tau);
    AbsAlgorithm<L,D>::x = x_initial;
    AbsAlgorithm<L,D>::primal_initial_set();
    AbsAlgorithm<L,D>::iter = 0;
    AbsAlgorithm<L,D>::restart_run_algorithm();
    /*
    AbsAlgorithm<L, D>::iter = 0;
    AbsAlgorithm<L, D>::time = 0.;
    D start;
    if(AbsAlgorithm<L, D>::result_iflog == true) {
        AbsAlgorithm<L, D>::result_log();
        Print("initial_left", AbsAlgorithm<L, D>::left);
        if(method_type == "restart"){
            while((AbsAlgorithm<L, D>::left > precision) && ( AbsAlgorithm<L,D>::iter < max_iter)) {
                start = clock();
                restart_step();
                AbsAlgorithm<L, D>::time += (clock() - start)/CLOCKS_PER_SEC;
                AbsAlgorithm<L, D>::result_log();
                Print("left", AbsAlgorithm<L, D>::left);
            }
        }else{
            sub_initial();
            while((AbsAlgorithm<L, D>::left > precision) && ( AbsAlgorithm<L,D>::iter < max_iter)) {
                start = clock();
                step();
                AbsAlgorithm<L, D>::time += (clock() - start)/CLOCKS_PER_SEC;
                AbsAlgorithm<L, D>::result_log();
                Print("left", AbsAlgorithm<L, D>::left);
            }
        }
    }else {
        if (method_type == "restart"){
            start = clock();
            while ((AbsAlgorithm<L, D>::left > precision) && (AbsAlgorithm<L, D>::iter < max_iter))
                restart_step();
            AbsAlgorithm<L, D>::time = (clock() - start) / CLOCKS_PER_SEC;
        } else {
            sub_initial();
            start = clock();
            while ((AbsAlgorithm<L, D>::left > precision) && (AbsAlgorithm<L, D>::iter < max_iter))
                step();
            AbsAlgorithm<L, D>::time = (clock() - start) / CLOCKS_PER_SEC;
        }
    }
    AbsAlgorithm<L, D>::logFile.close();
    */

}


template<typename L, typename D>
void APPROX<L, D>::result_print() {
    AbsAlgorithm<L, D>::loss->fun_value_compute(AbsAlgorithm<L, D>::fun_value, AbsAlgorithm<L, D>::data,
                                                AbsAlgorithm<L, D>::w, AbsAlgorithm<L, D>::x);
    AbsAlgorithm<L, D>::res_fun_value_compute(AbsAlgorithm<L, D>::fun_value);
    AbsAlgorithm<L, D>::loss->grad_compute(AbsAlgorithm<L, D>::grad, AbsAlgorithm<L, D>::data,
                                           AbsAlgorithm<L, D>::w, AbsAlgorithm<L, D>::x);
    AbsAlgorithm<L, D>::res_grad_compute(AbsAlgorithm<L, D>::grad);
    AbsAlgorithm<L, D>::grad_norm_update(AbsAlgorithm<L, D>::grad_norm);
    if(AbsAlgorithm<L, D>::method_type == "restart"){

        Print("ResAPPROX_K:", AbsAlgorithm<L, D>::K);
        Print("ResAPPROX_alpha:", AbsAlgorithm<L, D>::loss->alpha);
        //Print("ResAPPROX_b:", AbsAlgorithm<L, D>::loss->b);
        //Print("ResAPPROX_d:", AbsAlgorithm<L, D>::loss->d);
        Print("ResAPPROX_outer_iter:", AbsAlgorithm<L, D>::iter);
        Print("ResAPPROX_time:", AbsAlgorithm<L, D>::time);
        Print("ResAPPROX_initial_fun_value:",AbsAlgorithm<L, D>::initial_fun_value);
        Print("ResAPPROX_fun_value:", AbsAlgorithm<L, D>::fun_value);
        Print("ResAPPROX_initial_grad_norm:",AbsAlgorithm<L, D>::initial_grad_norm);
        Print("ResAPPROX_grad_norm:", AbsAlgorithm<L, D>::grad_norm);
        //Print("ResAPPROX_x:", AbsAlgorithm<L, D>::x);
        //Print("ResAPPROX_w:", AbsAlgorithm<L, D>::w);
        //Print("ResAPPROX_grad:", AbsAlgorithm<L, D>::grad);

    }else {
        Print("APPROX_outer_iter:", AbsAlgorithm<L, D>::iter);
        Print("APPROX_time:", AbsAlgorithm<L, D>::time);
        Print("APPROX_initial_fun_value:", AbsAlgorithm<L, D>::initial_fun_value);
        Print("APPROX_initial_grad_norm:", AbsAlgorithm<L, D>::initial_grad_norm);
        Print("APPROX_fun_value:", AbsAlgorithm<L, D>::fun_value);
        Print("APPROX_grad_norm:", AbsAlgorithm<L, D>::grad_norm);
        Print("APPROX_x:", AbsAlgorithm<L, D>::x);
        Print("APPROX_w:", AbsAlgorithm<L, D>::w);
        Print("APPROX_grad:", AbsAlgorithm<L, D>::grad);
    }
}



#endif //ALGORITHM_RESAPPROX_H
