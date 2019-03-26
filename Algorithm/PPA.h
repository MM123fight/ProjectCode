//
// Created by Lu Meng on 2019/2/20.
//

#ifndef ALGORITHM_PPA_H
#define ALGORITHM_PPA_H

#include "AbsAlgorithmSub.h"
#include "AbsAlgorithm.h"
#include "ParamPPA.h"
#include <stdio.h>



template<typename L, typename D>
class PPA:public AbsAlgorithmSub<L, D>{
private:
    L m,mi,n,tau;
    L sub_iter;
    L sub_max_iter = 1000;
    //primal_dual__gap,residual;
    D res_square, res_before;
    L t;
    D res_bound;
    D t_time;
    D outer_time;

    //parameters for PPA
    D sigma;
    D theta;
    D varepsilon_0, varepsilon;
    D delta_0,delta;

    D C,CK_param;
    D sub_precision,sub_alpha,right;
    D varepsilon_stage1;
    // varepsilon_stage2 = varepsilon_0; varepsilon_stage3 = varepsilon

    ParamPPA<D>* param;
    std::vector<D> zeros_n;

    bool stop;
    std::vector<D> sub_b;
    std::vector<D> sub_d;
    std::vector<D> x_before;
    std::vector<D> x_diff_over_sigma;
    std::vector<D> lambda_over_sigma;
    std::vector<D> lambda_over_sigma_before;
public:
    PPA(ProbData<L, D>* const data_inst, const unsigned int& sub_loss_type,
        const unsigned int& sub_reg_type,const unsigned int& sub_method_type,
        const std::string& data_file_path="");

    virtual ~PPA(){}
    void sub_step();
    void right_update();
    void sub_initial();
    void restart_sub_initial();

    void step();
    void restart_step();
    void solver(const std::vector<D>& x_initial,const std::vector<D>& lambda_initial,const L& sub_method_blocksize,
                const L& max_iter, const D& precision, ParamPPA<D>* const param,
                const bool& result_iflog = false, const std::string& method_type = "");
    //void res_initial();
    void left_res_initial();
    void left_res_compute();
    void result_print();
    void initial_check();
    void check();
    void restart_check();
    void initial_restart_check();
    void initial_t_check();
    void t_check();
    void initial_outer_check();
    void outer_check();


};

template<typename L, typename D>
PPA<L, D>::PPA(ProbData<L, D>* const data_inst, const unsigned int& sub_loss_type,
               const unsigned int& sub_reg_type,const unsigned int& sub_method_type,
               const std::string& data_file_path):
        AbsAlgorithmSub<L,D>(data_inst,sub_loss_type, sub_reg_type, sub_method_type,data_file_path){
    m = data_inst->m;
    mi = data_inst->mi;
    n = data_inst->n;
    zeros_n = std::vector<D>(n);
    sub_b = std::vector<D>(m);
    sub_d = std::vector<D>(n);
    x_before = std::vector<D>(n);
    x_diff_over_sigma = std::vector<D>(n);
    lambda_over_sigma = std::vector<D>(m);
    lambda_over_sigma_before = std::vector<D>(m);

    AbsAlgorithmSub<L, D>::x = std::vector<D>(n);
    AbsAlgorithmSub<L, D>::w = std::vector<D>(m);
    AbsAlgorithmSub<L, D>::lambda = std::vector<D>(m);
    AbsAlgorithmSub<L, D>::w_dual = std::vector<D>(n);

    D tmp = 2. * std::max<D>(l2norm_square<L,D>(data_inst->c),l2norm_square<L,D>(data_inst->b));
    std::vector<D> v(n);
    AbsAlgorithmSub<L, D>::set_loss();
    AbsAlgorithmSub<L, D>::loss->v_initial_set(v,data_inst,1);
    for (L row = 0; row < n; ++row)
        tmp += v[row];
    CK_param = sqrt(tmp);
    //Print("v",v);
    //data_inst->dataprint();
    //Print("CK_param",CK_param);
}


//fix sub_b, sub_d, varepsilon_k, delta_k,
template<typename L, typename D>
void PPA<L, D>::sub_initial() {
    sub_iter = 0;
}

template<typename L, typename D>
void PPA<L, D>::restart_sub_initial() {
    C = theta * CK_param * AbsAlgorithmSub<L, D>::left;
    //inv_update(varepsilon_stage1, param->varepsilon,param->inv_idx_varepsilon,AbsAlgorithmSub<L, D>::iter);
    varepsilon_stage1 *= param->exp_idx_varepsilon;
    AbsAlgorithmSub<L, D>::K = ceil(log(theta * CK_param)/param->beta);
    sub_alpha = 1./(sigma*sigma);
    t = 0;
    res_bound = AbsAlgorithmSub<L, D>::left;// res_bound = C*exp(-param->beta * K * t) for t = 0
    res_before = AbsAlgorithmSub<L, D>::left;
}

//PPA_step
template<typename L, typename D>
void PPA<L, D>::sub_step() {
    sub_precision = sub_alpha * varepsilon;
    //Print("sub_precision",sub_precision);
    for (L row = 0; row < m; ++row)
        sub_b[row] = AbsAlgorithmSub<L, D>::data->b[row] - lambda_over_sigma[row];
    for (L row = 0; row < n; ++row)
        sub_d[row] = AbsAlgorithmSub<L, D>::x[row] -sigma*AbsAlgorithmSub<L, D>::data->c[row];

    //initial_restart_check();

    AbsAlgorithmSub<L, D>::sub_method->set_loss(sub_alpha,sub_b,sub_d);
    AbsAlgorithmSub<L, D>::sub_method->solver(AbsAlgorithmSub<L, D>::x,tau,sub_max_iter,sub_precision,1,false,"restart");
    right_update();
    while(AbsAlgorithmSub<L, D>::sub_method->grad_norm > right) {
        AbsAlgorithmSub<L, D>::sub_method->restart_step();
        if(AbsAlgorithmSub<L, D>::sub_method->grad_norm <= sub_precision)
             right_update();
    }

    for (L row = 0; row < n; ++row)
        x_diff_over_sigma[row] = (AbsAlgorithmSub<L, D>::sub_method->x[row] - AbsAlgorithmSub<L, D>::x[row])/sigma;
    AbsAlgorithmSub<L, D>::x = AbsAlgorithmSub<L, D>::sub_method->x;
    D tmp;
    for (L row = 0; row < m; ++row){
        tmp = AbsAlgorithmSub<L, D>::sub_method->w[row];
        AbsAlgorithmSub<L, D>::w[row] = tmp - lambda_over_sigma[row];
        if((row < mi)&&(tmp < 0))
            tmp = 0;
        lambda_over_sigma[row] = tmp;
    }

    //restart_check();

    ++sub_iter;
}

template<typename L, typename D>
void PPA<L, D>::step() {
    //After every K steps, calculate the residual
    //The K should be a constant
    for (L i = 0; i < AbsAlgorithmSub<L, D>::K; ++i) {
        sigma *= param->exp_idx_theta;
        sub_alpha = 1./(sigma*sigma);
        inv_update<L,D>(varepsilon, varepsilon_0,param->inv_idx_varepsilon,sub_iter);
        inv_update<L,D>(delta,delta_0,param->inv_idx_delta,sub_iter);
        sub_step();
    }
    left_res_compute();
    //check();
    ++AbsAlgorithmSub<L, D>::iter;
}



template<typename L, typename D>
void PPA<L, D>::restart_step() {
    outer_time = 0.;
    D start = clock();
    restart_sub_initial();
    initial_outer_check();
    while( (AbsAlgorithmSub<L, D>::left <= res_bound)&&(AbsAlgorithmSub<L, D>::left <= res_before)
           &&(AbsAlgorithmSub<L, D>::left > AbsAlgorithmSub<L, D>::precision) ) {
        res_before = AbsAlgorithmSub<L, D>::left;
        x_before = AbsAlgorithmSub<L, D>::x;
        lambda_over_sigma_before = lambda_over_sigma;

        inv_update(varepsilon_0, varepsilon_stage1,param->inv_idx_varepsilon,t);
        //initial_t_check();
        //D start = clock();
        t_time = 0;
        //do PPA(K): with delta and sigma to be constant
        /********************************************************************************************/
        sub_iter = 0;
        while (sub_iter < AbsAlgorithmSub<L, D>::K) {
            inv_update<L,D>(varepsilon, varepsilon_0,param->inv_idx_varepsilon,sub_iter);
            sub_step();
        }
        left_res_compute();
        /********************************************************************************************/
        //t_time = (clock() - start)/(D)CLOCKS_PER_SEC;

        ++t;
        res_bound = C*exp(-param->beta * AbsAlgorithmSub<L, D>::K * t);
        //t_check();
    }

    if(AbsAlgorithmSub<L, D>::left > res_before) {
        AbsAlgorithmSub<L, D>::left = res_before;
        AbsAlgorithmSub<L, D>::x = x_before;
        lambda_over_sigma = lambda_over_sigma_before;
    }
    theta *= param->exp_idx_theta;
    sigma *= param->exp_idx_theta;
    ++AbsAlgorithmSub<L, D>::iter;
    outer_time = (clock() - start)/CLOCKS_PER_SEC;
    outer_check();
}



template<typename L, typename D>
void PPA<L, D>::solver(const std::vector<D>& x_initial, const std::vector<D>& lambda_initial, const L& sub_method_blocksize,
                       const L& max_iter, const D& precision,ParamPPA<D>* const param,
                       const bool& result_iflog, const std::string& method_type){
    tau = sub_method_blocksize;
    if(method_type == "restart")
        Print("PPA_method_type","restart PPA");
    else
        Print("PPA_method_type","PPA");
    if(result_iflog == true) {
        if(method_type == "restart") {
            AbsAlgorithmSub<L, D>::full_result_path = AbsAlgorithmSub<L, D>::data->data_file_path + "/result/ResPPA_subbs" +
                                                   to_string(tau);
        }else{
            AbsAlgorithmSub<L, D>::full_result_path = AbsAlgorithmSub<L, D>::data->data_file_path + "/result/PPA_subbs" +
                                                   to_string(tau);
        }
        AbsAlgorithmSub<L, D>::logfile_open();
    }

    AbsAlgorithmSub<L, D>::x = x_initial;
    AbsAlgorithmSub<L, D>::lambda = lambda_initial;
    AbsAlgorithmSub<L, D>::max_iter = max_iter;
    AbsAlgorithmSub<L, D>::precision = precision;
    this->param = param;
    AbsAlgorithmSub<L, D>::result_iflog = result_iflog;
    AbsAlgorithmSub<L, D>::method_type = method_type;

    AbsAlgorithmSub<L, D>::K = 1;
    varepsilon_0 = param->varepsilon;
    varepsilon_stage1 = param->varepsilon;
    delta_0 = param->delta;
    delta = param->delta;
    theta = param->theta;
    sigma = param->theta*param->alpha;
    for (L row = 0; row < m; ++row) {
        lambda_over_sigma[row] = lambda_initial[row]/sigma;
    }
    //w = Ax-b w_dual = AT*lambda+c
    AbsAlgorithmSub<L, D>::loss->w_initial(AbsAlgorithmSub<L, D>::w,AbsAlgorithmSub<L, D>::data, x_initial,
                                          AbsAlgorithmSub<L, D>::data->b);
    AbsAlgorithmSub<L, D>::loss->w_dual_initial(AbsAlgorithmSub<L, D>::w_dual,AbsAlgorithmSub<L, D>::data,
                                                 lambda_initial,AbsAlgorithmSub<L, D>::data->c);
    left_res_initial();
    initial_check();
    AbsAlgorithmSub<L, D>::restart_run_algorithm();

}

template<typename L, typename D>
void PPA<L, D>::left_res_initial(){
    AbsAlgorithmSub<L, D>::sub_method->grad = zeros_n;
    x_diff_over_sigma = zeros_n;
    left_res_compute();
    D tmp;
    for (L row = 0; row < n; ++row) {
        tmp = AbsAlgorithmSub<L, D>::w_dual[row];
        res_square += tmp * tmp;
    }
    AbsAlgorithmSub<L, D>::left = sqrt(res_square);
}


template<typename L, typename D>
void PPA<L, D>::left_res_compute() {
    //primal_dual_gap
    D tmp;
    AbsAlgorithmSub<L, D>::fun_value = 0;
    for (L row = 0; row < n; ++row)
        AbsAlgorithmSub<L, D>::fun_value += AbsAlgorithmSub<L, D>::data->c[row] * AbsAlgorithmSub<L, D>::x[row];

    AbsAlgorithmSub<L, D>::primal_dual_gap = 0;
    for (L row = 0; row < m; ++row)
        AbsAlgorithmSub<L, D>::primal_dual_gap += AbsAlgorithmSub<L, D>::data->b[row]*lambda_over_sigma[row];
    AbsAlgorithmSub<L, D>::primal_dual_gap *= sigma;
    AbsAlgorithmSub<L, D>::primal_dual_gap += AbsAlgorithmSub<L, D>::fun_value;
    tmp = AbsAlgorithmSub<L, D>::primal_dual_gap;
    res_square = tmp*tmp;
    //primal_res
    for (L row = 0; row < m; ++row){
        tmp = AbsAlgorithmSub<L, D>::w[row];
        if((row >= mi)||(tmp > 0))
            res_square += tmp*tmp;
    }
    //dual_res
    for (L row = 0; row < n; ++row) {
        tmp = sigma * AbsAlgorithmSub<L, D>::sub_method->grad[row] - x_diff_over_sigma[row];
        res_square += tmp * tmp;
    }
    AbsAlgorithmSub<L, D>::left = sqrt(res_square);
}

template<typename L, typename D>
void PPA<L, D>::right_update() {
    right = 0;
    D tmp;
    for (L row = 0; row < n; ++row) {
        tmp = (AbsAlgorithmSub<L, D>::sub_method->x[row] - AbsAlgorithmSub<L, D>::x[row])/sigma;
        right += tmp*tmp;
    }
    //AbsAlgorithmSub<L, D>::sub_method->result_print();

    for (L row = 0; row < m; ++row) {
        tmp = AbsAlgorithmSub<L, D>::sub_method->w[row];
        if((row < mi) && (tmp < 0)) {
            tmp = 0;
        }
        tmp -= lambda_over_sigma[row];
        right += tmp * tmp;
    }
    right = sqrt(right);
    right *= (delta/sigma);
    /*
    Print("sub_outer_iter",AbsAlgorithmSub<L, D>::sub_method->iter);
    Print("sub_grad_norm",AbsAlgorithmSub<L, D>::sub_method->grad_norm);
    Print("right", right);
     */
}

template<typename L, typename D>
void PPA<L, D>::result_print() {
    Print();
    if(AbsAlgorithmSub<L, D>::method_type == "restart"){
        Print("ResPPA_outer_iter:", AbsAlgorithmSub<L, D>::iter);
        Print("ResPPA_time:", AbsAlgorithmSub<L, D>::time);
        Print("ResPPA_obj_fun_value:", AbsAlgorithmSub<L, D>::fun_value);
        Print("ResPPA_primal_dual_gap",AbsAlgorithmSub<L, D>::primal_dual_gap);
        Print("ResPPA_residual:", AbsAlgorithmSub<L, D>::left);
        Print("ResPPA_x:", AbsAlgorithmSub<L, D>::x);
        Print("ResPPA_w:", AbsAlgorithmSub<L, D>::w);
        Print("ResPPA_lambda:", AbsAlgorithmSub<L, D>::lambda);
        Print("ResPPA_primal_dual_gap:",AbsAlgorithmSub<L, D>::primal_dual_gap);

    }else {
        Print("PPA_outer_iter:", AbsAlgorithmSub<L, D>::iter);
        Print("PPA_time:", AbsAlgorithmSub<L, D>::time);
        Print("PPA_obj_fun_value:", AbsAlgorithmSub<L, D>::fun_value);
        Print("PPA_primal_dual_gap",AbsAlgorithmSub<L, D>::primal_dual_gap);
        Print("PPA_residual:", AbsAlgorithmSub<L, D>::left);
        Print("PPA_x:", AbsAlgorithmSub<L, D>::x);
        Print("PPA_w:", AbsAlgorithmSub<L, D>::w);
        Print("PPA_lambda:", AbsAlgorithmSub<L, D>::lambda);
        Print("PPA_primal_dual_gap:",AbsAlgorithmSub<L, D>::primal_dual_gap);
        Print("PPA_residual:",sqrt(res_square));
    }
}

template<typename L, typename D>
void PPA<L, D>::initial_check() {
    Print();
    Print("initial_check start!");
    Print("x:", AbsAlgorithmSub<L, D>::x);
    Print("w:", AbsAlgorithmSub<L, D>::w);
    Print("lambda:", lambda_over_sigma);
    Print("Obj_fun_value:", AbsAlgorithmSub<L, D>::fun_value);
    Print("primal_dual_gap",AbsAlgorithmSub<L, D>::primal_dual_gap);
    Print("residual",AbsAlgorithmSub<L, D>::left);
    Print("initial_check end!");
    Print();

}

template<typename L, typename D>
void PPA<L, D>::check() {
    Print("check start!");
    //For the algorithm, we don't need the value of w_dual, thus, we do not calculate.
    Print("sub_alpha", sub_alpha);
    Print("sub_b", sub_b);
    Print("sub_d", sub_d);


    Print("sub_x",AbsAlgorithmSub<L, D>::sub_method->x);
    Print("x_diff_over_sigma",x_diff_over_sigma);
    Print("sub_w",AbsAlgorithmSub<L, D>::sub_method->w);
    Print("sub_grad", AbsAlgorithmSub<L, D>::sub_method->grad);
    Print("sub_left", AbsAlgorithmSub<L, D>::sub_method->left);
    Print("sub_right",std::min<D>(right,sub_precision));

    Print("x:", AbsAlgorithmSub<L, D>::x);
    Print("w:", AbsAlgorithmSub<L, D>::w);
    Print("sigma:", sigma);
    Print("lambda_over_sigma:", lambda_over_sigma);
    Print("Obj_fun_value:", AbsAlgorithmSub<L, D>::fun_value);
    Print("primal_dual_gap",AbsAlgorithmSub<L, D>::primal_dual_gap);
    Print("residual",AbsAlgorithmSub<L, D>::left);

    Print("check end!");
    Print();
}

template <typename L, typename D>
void PPA<L, D>::initial_restart_check(){
    Print("--------------------------------------------------------");
    Print("sub_iter", sub_iter);
    Print("sub_alpha", sub_alpha);
    Print("varepsilon_stage3", varepsilon);

    std::string sub_b_path = AbsAlgorithm<L, D>::data->data_file_path + "/result/sub_b";
    VectorWrite(sub_b_path,sub_b);
    std::string sub_d_path = AbsAlgorithm<L, D>::data->data_file_path + "/result/sub_d";
    VectorWrite(sub_d_path,sub_d);

    std::string initial_x_path = AbsAlgorithm<L, D>::data->data_file_path + "/result/initial_x";
    std::string initial_lambda_over_sigma_path = AbsAlgorithm<L, D>::data->data_file_path + "/result/initial_lambda_over_sigma";
    VectorWrite(initial_x_path,AbsAlgorithmSub<L, D>::x);
    VectorWrite(initial_lambda_over_sigma_path,lambda_over_sigma);
    //Print("sub_b",sub_b);
    //Print("sub_d",sub_d);
    //Print("initial_x", AbsAlgorithmSub<L, D>::x);
    //Print("initial_lambda_over_sigma", lambda_over_sigma);
}

template <typename L, typename D>
void PPA<L, D>::restart_check() {
    //For the algorithm, we don't need the value of w_dual, thus, we do not calculate.
    AbsAlgorithmSub<L, D>::sub_method->result_print();
    Print("right", right);
    std::string x_path = AbsAlgorithm<L, D>::data->data_file_path + "/result/x";
    std::string lambda_over_sigma_path = AbsAlgorithm<L, D>::data->data_file_path + "/result/lambda_over_sigma";
    VectorWrite(x_path,AbsAlgorithmSub<L, D>::x);
    VectorWrite(lambda_over_sigma_path,lambda_over_sigma);
    //Print("x",AbsAlgorithmSub<L, D>::x);
    //Print("lambda_over_sigma", lambda_over_sigma);
    Print("--------------------------------------------------------");

}

template <typename L, typename D>
void PPA<L, D>::initial_t_check(){
    Print("--------------------------------------------------------");
    Print("t", t);
    Print("varepsilon_stage2", varepsilon_0);
    Print("initial_obj_fun_value", AbsAlgorithmSub<L, D>::fun_value);
    Print("initial_left",AbsAlgorithmSub<L, D>::left);
    std::string initial_x_path = AbsAlgorithm<L, D>::data->data_file_path + "/result/initial_x";
    std::string initial_lambda_over_sigma_path = AbsAlgorithm<L, D>::data->data_file_path + "/result/initial_lambda_over_sigma";
    VectorWrite(initial_x_path,AbsAlgorithmSub<L, D>::x);
    VectorWrite(initial_lambda_over_sigma_path,lambda_over_sigma);
}


template <typename L, typename D>
void PPA<L, D>::t_check() {
    Print("obj_fun_value", AbsAlgorithmSub<L, D>::fun_value);
    Print("left",AbsAlgorithmSub<L, D>::left);
    Print("res_bound", res_bound);
    Print("t_time",t_time);
    std::string x_path = AbsAlgorithm<L, D>::data->data_file_path + "/result/x";
    std::string lambda_over_sigma_path = AbsAlgorithm<L, D>::data->data_file_path + "/result/lambda_over_sigma";
    VectorWrite(x_path,AbsAlgorithmSub<L, D>::x);
    VectorWrite(lambda_over_sigma_path,lambda_over_sigma);
    //Print("x",AbsAlgorithmSub<L, D>::x);
    //Print("lambda_over_sigma", lambda_over_sigma);
    Print("--------------------------------------------------------");

}


template <typename L, typename D>
void PPA<L, D>::initial_outer_check() {
    printf("\033[31m---------------------------------------------------------------\033[0m \n");
    Print("outer_iter", AbsAlgorithmSub<L, D>::iter);
    Print("C", C);
    Print("K", AbsAlgorithmSub<L, D>::K);
    Print("sigma", sigma);
    Print("varepsilon_stage1", varepsilon_stage1);
}

template <typename L, typename D>
void PPA<L, D>::outer_check() {
    Print("left", AbsAlgorithmSub<L, D>::left);
    Print("fun_value",AbsAlgorithmSub<L, D>::fun_value);
    Print("outer_time", outer_time);
    printf("\033[31m---------------------------------------------------------------\033[0m \n");
    Print();
}


#endif //ALGORITHM_PPA_H
