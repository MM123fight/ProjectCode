#ifndef ALGORITHM_PCDM_H
#define ALGORITHM_PCDM_H

#include "AbsAlgorithm.h"
#include "AbsAlgorithm.h"

template<typename L, typename D>
class PCDM:public AbsAlgorithm<L, D>{
private:
    L m,n,tau;
    std::vector<D> delta_x;
    std::vector<L> cord;
    std::vector<D> cord_grad;
    D theta;
    D precision;
    gsl_rng *gsl_rng_r;
public:
    PCDM(ProbData<L, D>* const data_inst, const unsigned int& loss_type_value, const unsigned int& reg_type_value,
         const std::string& data_file_path_value="");
    ~PCDM() {
        gsl_rng_free(gsl_rng_r);
    }
    void sub_step();
    void step();
    void solver(const std::vector<D>& x_initial, const L& blocksize, const L& max_iter,
                const D& precision_value, const unsigned int& stop_type = 1, const bool& result_iflog_value = false,
                const std::string& method_type = "");
    void initial_print();
    void result_print();


};
template<typename L, typename D>
PCDM<L, D>::PCDM(ProbData<L, D> *const data_inst, const unsigned int &loss_type_value,const unsigned int &reg_type_value,
                 const std::string& data_file_path_value):
        AbsAlgorithm<L, D>(data_inst, loss_type_value,reg_type_value,data_file_path_value) {
    gsl_rng_env_setup();
    gsl_rng_r = gsl_rng_alloc(gsl_rng_default);
    m = data_inst->m;
    n = data_inst->n;
    AbsAlgorithm<L,D>::grad = std::vector<D>(n);
    AbsAlgorithm<L,D>::w = std::vector<D>(m);
    AbsAlgorithm<L,D>::x = std::vector<D>(n);
}


template<typename L, typename D>
void PCDM<L, D>::sub_step() {
    AbsAlgorithm<L,D>::loss->cord_grad_update(cord_grad,AbsAlgorithm<L,D>::data,AbsAlgorithm<L,D>::w,AbsAlgorithm<L,D>::x,tau,cord);
    AbsAlgorithm<L,D>::reg_cord_cal(AbsAlgorithm<L,D>::x,delta_x,1.,cord_grad,tau,cord);
    AbsAlgorithm<L,D>::loss->w_update(AbsAlgorithm<L,D>::w,AbsAlgorithm<L,D>::data,delta_x,tau,cord);
}

template<typename L, typename D>
void PCDM<L, D>::step() {

    D start = clock();
    for (L i = 0; i < AbsAlgorithm<L, D>::K; ++i) {
        SampleGen(cord,tau,n,gsl_rng_r);
        sub_step();
    }
    AbsAlgorithm<L, D>::left_value_compute();
    ++AbsAlgorithm<L,D>::iter;
    AbsAlgorithm<L, D>::time += (D)(clock() - start)/CLOCKS_PER_SEC;
    if(AbsAlgorithm<L, D>::result_iflog == true)
        AbsAlgorithm<L, D>::result_log();
    Print("left", AbsAlgorithm<L, D>::left);
}


template<typename L, typename D>
void PCDM<L, D>::solver(const std::vector<D>& x_initial, const L& blocksize, const L& max_iter,
                        const D& precision_value, const unsigned int& stop_type_value, const bool& result_iflog_value,
                        const std::string& method_type){
    AbsAlgorithm<L, D>::check_blocksize(blocksize);
    tau = blocksize;
    AbsAlgorithm<L, D>::result_iflog = result_iflog_value;
    if(AbsAlgorithm<L, D>::result_iflog == true) {
        AbsAlgorithm<L, D>::full_result_path = AbsAlgorithm<L, D>::data->data_file_path + "/result/PCDM_bs" + to_string(tau)
                                       + "_alpha" + to_string(AbsAlgorithm<L, D>::loss->alpha);
        AbsAlgorithm<L, D>::logfile_open();
    }

    AbsAlgorithm<L, D>::stop_type = stop_type_value;
    precision = precision_value;
    AbsAlgorithm<L, D>::set_v_mu(tau);
    AbsAlgorithm<L, D>::Krcd_compute(AbsAlgorithm<L, D>::K,tau);

    cord = std::vector<L>(tau);
    cord_grad = std::vector<D>(tau);
    delta_x = std::vector<D>(tau);

    //Print("K:",K);

    AbsAlgorithm<L,D>::x.assign(x_initial.begin(), x_initial.end());
    AbsAlgorithm<L,D>::primal_initial_set();
    AbsAlgorithm<L,D>::iter = 0;
    AbsAlgorithm<L,D>::time = 0;
    if(AbsAlgorithm<L, D>::result_iflog == true)
        AbsAlgorithm<L, D>::result_log();
    Print("intial_left", AbsAlgorithm<L,D>::left);
    while((AbsAlgorithm<L,D>::left > precision) && ( AbsAlgorithm<L,D>::iter < max_iter)){
        step();
    }
    AbsAlgorithm<L, D>::logFile.close();

}


template<typename L, typename D>
void PCDM<L, D>::initial_print() {
    Print("PCDM_initial_fun_value:",AbsAlgorithm<L, D>::initial_fun_value);
    Print("PCDM_initial_grad_norm:",AbsAlgorithm<L, D>::initial_grad_norm);
}

template<typename L, typename D>
void PCDM<L, D>::result_print() {
    Print("PCDM_outer_iter:", AbsAlgorithm<L, D>::iter);
    Print("PCDM_time:", AbsAlgorithm<L, D>::time);
    AbsAlgorithm<L,D>::fun_value_compute(AbsAlgorithm<L,D>::fun_value, AbsAlgorithm<L,D>::w, AbsAlgorithm<L,D>::x);
    Print("PCDM_fun_value:",AbsAlgorithm<L, D>::fun_value);
    AbsAlgorithm<L,D>::grad_compute(AbsAlgorithm<L,D>::grad,AbsAlgorithm<L,D>::w, AbsAlgorithm<L,D>::x);
    AbsAlgorithm<L,D>::grad_norm_update(AbsAlgorithm<L,D>::grad_norm);
    Print("PCDM_grad_norm:",AbsAlgorithm<L, D>::grad_norm);
    Print("PCDM_x:",AbsAlgorithm<L, D>::x);
    Print("PCDM_w:",AbsAlgorithm<L, D>::w);
    Print("PCDM_grad:",AbsAlgorithm<L, D>::grad);
}



#endif //ALGORITHM_PCDM_H
