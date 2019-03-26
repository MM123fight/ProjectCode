#ifndef ALGORITHM_NU_ACDM_H
#define ALGORITHM_NU_ACDM_H

#include "AbsAlgorithm.h"
#include "AbsAlgorithm.h"
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_cblas.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>


template<typename L, typename D>
class NU_ACDM:public AbsAlgorithm<L, D>{
private:
    D beta = 0;
    D gamma, S_gamma, tau, eta;
    L m, n, cord;
    D grad_cord, sigma;
    D eta_1, eta_2, eta_3;
    D tmp, tmp1, tmp2, tmp_value;
    L tmp_idx;
    std::vector<D> M_base;
    std::vector<D> M;
    std::vector<D> M_inv_base;
    std::vector<D> M_inv;
    std::vector<D> Matrix_tmp;

    //y_bar = M*y; z_bar = M*z;
    //w_y_bar = M*w_y; w_z_bar = M*w_z;
    std::vector<D> y_bar;
    std::vector<D> z_bar;
    std::vector<D> w_y_bar;
    std::vector<D> w_z_bar;
    std::vector<D> Lip_gamma;
    std::vector<D> Lip_beta;
    std::vector<D> p;
    std::vector<D> p_sum;
    D eta_cord_y;
    D eta_cord_z;
    gsl_rng *gsl_rng_r;

public:
    NU_ACDM(ProbData<L, D>* const data_inst, const unsigned int& loss_type, const unsigned int& reg_type,
            const std::string& data_file_path="");
    ~NU_ACDM() {
        gsl_rng_free(gsl_rng_r);
    }
    void sub_initial();
    void sub_step();
    void step();
    void solver(const std::vector<D>& x_initial, const L& max_iter, const D& precision,
                const unsigned int& stop_type, const bool& result_iflog);
    void left_value_compute();
    void initial_print();
    void result_print();


};

template<typename L, typename D>
NU_ACDM<L, D>::NU_ACDM(ProbData<L, D>* const data_inst, const unsigned int& loss_type, const unsigned int& reg_type,
                       const std::string& data_file_path):
        AbsAlgorithm<L, D>(data_inst, loss_type, reg_type, data_file_path) {
    gsl_rng_env_setup();
    gsl_rng_r = gsl_rng_alloc(gsl_rng_default);
    m = data_inst->m;
    n = data_inst->n;
    M_base = std::vector<D>(4);
    M = std::vector<D>(4);
    M_inv_base = std::vector<D>(4);
    M_inv = std::vector<D>(4);
    Matrix_tmp = std::vector<D>(4);

    y_bar = std::vector<D>(n);
    z_bar = std::vector<D>(n);
    w_y_bar = std::vector<D>(m);
    w_z_bar = std::vector<D>(m);
    p = std::vector<D>(n);
    p_sum = std::vector<D>(n);
    Lip_gamma = std::vector<D>(n);
    Lip_beta = std::vector<D>(n);
    AbsAlgorithm<L,D>::Lip = std::vector<D>(n);
    AbsAlgorithm<L,D>::grad = std::vector<D>(n);
    AbsAlgorithm<L,D>::w = std::vector<D>(m);
    AbsAlgorithm<L,D>::x = std::vector<D>(n);

}

template<typename L, typename D>
void NU_ACDM<L, D>::sub_initial() {
    sigma = AbsAlgorithm<L, D>::loss->alpha;
    gamma = 0.5*(1-beta);
    AbsAlgorithm<L,D>::loss->Lip_set(AbsAlgorithm<L,D>::Lip,AbsAlgorithm<L,D>::data);
    Print("Lip",AbsAlgorithm<L,D>::Lip );
    for (L row = 0; row < n; ++row) {
        Lip_gamma[row] = pow(AbsAlgorithm<L, D>::Lip[row], gamma);
        Lip_beta[row] = pow(AbsAlgorithm<L, D>::Lip[row], beta);
    }

    S_gamma = 0.;
    for (L row = 0; row < n; ++row)
        S_gamma += Lip_gamma[row];

    for (L row = 0; row < n; ++row)
        p[row] = Lip_gamma[row]/S_gamma;
    for (L row = 0; row < n-1; ++row) {
        p_sum[row + 1] = p_sum[row] + p[row];
    }
    Print("CheckLip_2",AbsAlgorithm<L,D>::Lip );
    tau = 2. / (1 + sqrt(4 * S_gamma * S_gamma/sigma + 1));
    eta = 1. / (tau * S_gamma*S_gamma);
    eta_3 = eta/(1+eta*sigma);
    eta_1 = tau * eta_3;
    eta_2 = 1 - tau;
    Print("eta", eta);
    Print("eta_1", eta_1);
    Print("eta_2", eta_2);
    Print("eta_3", eta_3);

    M_base[1] = tau/(1+eta*sigma);
    M_base[0] = 1- M_base[1];
    M_base[3] = 1./(1+eta*sigma);
    M_base[2] = 1-M_base[3];
    M[0] = 1; M_inv[0] = 1;
    M[1] = 0; M_inv[1] = 0;
    M[2] = 0; M_inv[2] = 0;
    M[3] = 1; M_inv[3] = 1;
    Print("M", M_base);

    D determint = M_base[0] * M_base[3] - M_base[1] * M_base[2];
    Print("determint", determint);
    M_inv_base[0] = M_base[3]/determint;
    M_inv_base[1] = - M_base[1]/determint;
    M_inv_base[2] = - M_base[2]/determint;
    M_inv_base[3] = M_base[0]/determint;
    Print("M_inv_base",M_inv_base);


    //y = x; z = x; w_y = A*y; w_z = A*z;
    y_bar = AbsAlgorithm<L, D>::x;
    z_bar = AbsAlgorithm<L, D>::x;
    AbsAlgorithm<L, D>::loss->w_initial(w_y_bar, AbsAlgorithm<L, D>::data, AbsAlgorithm<L, D>::x);
    w_z_bar = w_y_bar;
    Print("check_Lip",AbsAlgorithm<L,D>::Lip );

}

template<typename L, typename D>
void NU_ACDM<L, D>::sub_step() {

    AbsAlgorithm<L,D>::loss->cord_grad_update(grad_cord,AbsAlgorithm<L,D>::data, M[0], M[1], w_y_bar, y_bar, w_z_bar, z_bar, cord);
    Print("cord", cord);
    Print("grad_cord", grad_cord);
    tmp = p[cord];
    Print("Lip", AbsAlgorithm<L, D>::Lip);
    tmp1 = (eta_1/tmp+ eta_2/AbsAlgorithm<L, D>::Lip[cord])*grad_cord;
    tmp2 = (eta_3/tmp)*grad_cord;
    Print("tmp",tmp);
    Print("tmp1",tmp1);
    Print("tmp2",tmp2);

    MtimesM_s2(Matrix_tmp,M,M_base);
    M = Matrix_tmp;
    Print("M", M);
    MtimesM_s2(Matrix_tmp,M_inv,M_inv_base);
    M_inv = Matrix_tmp;
    Print("M_inv", M_inv);

    eta_cord_y = M_inv[0]*tmp1 + M_inv[1]*tmp2;
    eta_cord_z = M_inv[2]*tmp1 + M_inv[3]*tmp2;
    Print("eta_cord_y", eta_cord_y);
    Print("eta_cord_z", eta_cord_z);

    y_bar[cord] -= eta_cord_y;
    z_bar[cord] -= eta_cord_z;
    Print("y_bar", y_bar);
    Print("z_bar", z_bar);
    std::vector<D> y(n);
    std::vector<D> z(n);
    for (int i = 0; i < n; ++i) {
        y[i] = M[0] * y_bar[i] + M[1] * z_bar[i];
        z[i] = M[2] * y_bar[i] + M[3] * z_bar[i];
    }
    Print("y",y);
    Print("z",z);
    for (L row = AbsAlgorithm<L, D>::data->AT_row_ptr[cord]; row < AbsAlgorithm<L, D>::data->AT_row_ptr[cord+1]; ++row) {
        tmp_idx =  AbsAlgorithm<L, D>::data->AT_col_idx[row];
        tmp_value =  AbsAlgorithm<L, D>::data->AT_value[row];
        w_y_bar[tmp_idx] -= eta_cord_y * tmp_value;
        w_z_bar[tmp_idx] -= eta_cord_z * tmp_value;
    }
    // M update M = M * M_base;
}

template<typename L, typename D>
void NU_ACDM<L, D>::step() {
    for (L i = 0; i < AbsAlgorithm<L, D>::K; ++i) {
        cord = RandSample(n,p_sum,gsl_rng_r) - 1;
        sub_step();
    }
    left_value_compute();
    //Print("grad_norm:",AbsAlgorithm<L, D>::grad_norm);
    ++AbsAlgorithm<L,D>::iter;
    Print(AbsAlgorithm<L,D>::iter);
}


template<typename L, typename D>
void NU_ACDM<L, D>::solver(const std::vector<D>& x_initial, const L& max_iter,
                           const D& precision, const unsigned int& stop_type, const bool& result_iflog){
    AbsAlgorithm<L, D>::result_iflog = result_iflog;
    if(AbsAlgorithm<L, D>::result_iflog == true) {
        AbsAlgorithm<L, D>::full_result_path = AbsAlgorithm<L, D>::data->data_file_path + "/result/NU_ACDM_alpha" + to_string(AbsAlgorithm<L, D>::loss->alpha);
        AbsAlgorithm<L, D>::logfile_open();
    }
    AbsAlgorithm<L, D>::stop_type = stop_type;
    AbsAlgorithm<L, D>::precision = precision;
    AbsAlgorithm<L, D>::max_iter = max_iter;
    AbsAlgorithm<L, D>::set_v_mu(1);
    AbsAlgorithm<L, D>::Krcd_compute(AbsAlgorithm<L, D>::K,1);

    AbsAlgorithm<L, D>::x = x_initial;
    AbsAlgorithm<L,D>::primal_initial_set();
    sub_initial();
    AbsAlgorithm<L,D>::run_algorithm();
}

template<typename L, typename D>
void NU_ACDM<L, D>::left_value_compute() {
    Print("check");
    switch (AbsAlgorithm<L, D>::stop_type) {
        case 1:
            AbsAlgorithm<L,D>::loss->grad_compute(AbsAlgorithm<L,D>::grad,AbsAlgorithm<L,D>::data,
                                                  tau, M, w_y_bar, y_bar, w_z_bar, z_bar);
            AbsAlgorithm<L,D>::res_grad_compute(AbsAlgorithm<L,D>::grad);
            AbsAlgorithm<L,D>::grad_norm_update(AbsAlgorithm<L,D>::grad_norm);
            AbsAlgorithm<L,D>::left = AbsAlgorithm<L,D>::grad_norm;
            break;
        default:
            Print("There is no such stop type!");
            break;
    }
}


template<typename L, typename D>
void NU_ACDM<L, D>::initial_print() {
    Print("NU_ACDM_initial_fun_value:",AbsAlgorithm<L, D>::initial_fun_value);
    Print("NU_ACDM_initial_grad_norm:",AbsAlgorithm<L, D>::initial_grad_norm);
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
