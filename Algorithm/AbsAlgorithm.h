//Lu Meng
#ifndef ALGORITHM_ABSALGORITHM_H
#define ALGORITHM_ABSALGORITHM_H

#include "../Problem/ProblemHeader.h"
#include "../LossFun/LossFunHeader.h"
#include "../use/useheader.h"
#include <algorithm>
#include <cmath>

template<typename L, typename D>
class AbsAlgorithm{
protected:
    ProbData<L, D>* data = NULL;
    AbsLoss<L, D>* loss = NULL;
    //random algorithm
    D mu;
    std::vector<D> v;
    std::vector<D> Lip;
    unsigned int loss_type;
    unsigned int stop_type;
    unsigned int reg_type;
    std::ofstream logFile;
    std::string data_file_path;
    std::string full_result_path;
    D left;
    bool result_iflog;
    std::vector<L> full_block;
    D reg_alpha = 1.;

public:
    L K;
    D initial_grad_norm,grad_norm;
    D initial_fun_value, fun_value;
    D initial_residual,residual;
    D initial_primal_dual_grap,primal_dual_gap;
    std::vector<D> grad;
    std::vector<D> w;//res_primal;
    std::vector<D> x;
    std::vector<D> lambda;
    std::vector<D> res_dual;
    L iter;
    D time;
    D fun_value_opt = 0;
    AbsAlgorithm(ProbData<L, D>* const data_inst, const unsigned int& loss_type_value, const unsigned int& reg_type_value,
                 const std::string& data_file_path_value=""):
            data(data_inst),loss_type(loss_type_value), reg_type(reg_type_value), data_file_path(data_file_path_value){
        L n = data->n;
        full_block = std::vector<L>(n);
        for (L row = 0; row < n; ++row)
            full_block[row] = row;
    }

    virtual ~AbsAlgorithm() {
        delete loss;
        if(logFile.is_open()) {
            logFile.close();
        }
    }
    void set_v_mu(const L& blocksize);
    void set_loss();
    void set_loss(const std::vector<D>&b);
    void set_loss(const D& alpha, const std::vector<D>&b, const std::vector<D>& d);
    void set_reg(const D& alpha);
    void grad_norm_update(D & grad_norm);
    void before_test(const D& value, const D& value_before, std::vector<D>& vec1, std::vector<D>& vec1_before,
                     std::vector<D>& vec2, std::vector<D>& vec2_before);
    void check_blocksize(const L& blocksize);
    void left_value_compute();
    void logfile_open();
    void result_log();
    void Krcd_compute(L& K, const L& blocksize);
    void primal_initial_set();
    void reg_cord_cal(std::vector<D> &x,std::vector<D>& delta_x, const D& coefficient ,
                      const std::vector<D> &cord_grad, const L& blocksize, const std::vector<L>& cord) ;
    virtual void sub_initial(){}
    virtual void sub_step(){}
    virtual void step(){}
    virtual void restart_step(){}
    virtual void solver(const std::vector<D>& x_initial, const L& blocksize, const L& max_iter,
                const D& precision_value, const unsigned int& stop_type = 1, const bool& result_iflog_value = false,
                const std::string& method_type = ""){}
    virtual void right_update(){}
    virtual void res_update(){}
    virtual void initial_print(){}
    virtual void result_print(){}

    void grad_compute(std::vector<D>& grad,const std::vector<D>& w, const std::vector<D>& x);

    void grad_sum_compute(std::vector<D>& grad_sum,const std::vector<D>& w_y, const std::vector<D>& y,
                          const std::vector<D>& w_z, const std::vector<D>& z, const D& theta_y);
    void fun_value_compute(D& fun_value,const std::vector<D>& w, const std::vector<D>& x);


};


template<typename L, typename D>
void AbsAlgorithm<L, D>::set_v_mu(const L& blocksize){
    loss->v_set(v,data,blocksize);
    //nested dependent names need "typename" to verify.
    typename std::vector<D>::iterator maxPosition = std::max_element(v.begin(),v.end());
    mu = loss->alpha / (*maxPosition);
}

template<typename L, typename D>
void AbsAlgorithm<L, D>::set_loss() {
    loss = new AbsLoss<L, D>();
}

template<typename L, typename D>
void AbsAlgorithm<L, D>::set_loss(const std::vector<D>&b) {
    if(loss == NULL)
        loss = new squared_loss<L, D>(b);
    else{
        loss->set_b(b);
    }
}

template<typename L, typename D>
void AbsAlgorithm<L, D>::set_loss(const D& alpha, const std::vector<D>&b, const std::vector<D>& d) {
    if(loss == NULL)
        loss = new plus_squared_loss<L, D>(alpha, b, d);
    else{
        loss->set_b(b);
        loss->set_d(d);
        loss->alpha = alpha;
    }
}

template<typename L, typename D>
void AbsAlgorithm<L, D>::set_reg(const D& alpha) {
    reg_alpha = alpha;
}


template<typename L, typename D>
void AbsAlgorithm<L, D>::grad_norm_update(D & grad_norm) {
    D tmp;
    grad_norm = 0;
    for (L row = 0; row < data->n; ++row) {
        tmp = grad[row];
        grad_norm +=  tmp*tmp;
    }
    grad_norm = sqrt(grad_norm);
}

template<typename L, typename D>
void AbsAlgorithm<L, D>::before_test(const D& value, const D& value_before, std::vector<D>& vec1, std::vector<D>& vec1_before,
                                     std::vector<D>& vec2, std::vector<D>& vec2_before){

    if( value > value_before){
        for (L row = 0; row < vec1.size(); ++row)
            vec1[row] = vec1_before[row];
        for (L row = 0; row < vec2.size(); ++row)
            vec2[row] = vec2_before[row];
    }else{
        for (L row = 0; row < vec1.size(); ++row)
            vec1_before[row] =  vec1[row];
        for (L row = 0; row < vec2.size(); ++row)
            vec2_before[row] =  vec2[row];
    }
}

template<typename L, typename D>
void AbsAlgorithm<L, D>::left_value_compute() {
    switch (stop_type) {
        case 1:
            grad_compute(grad, w, x);
            grad_norm_update(grad_norm);
            left = grad_norm;
            break;
        case 2:
            loss->fun_value_compute(fun_value, data, w, x);
            left = fun_value - fun_value_opt;
            break;
        default:
            Print("There is no such stop type!");
            break;
    }
}
template<typename L, typename D>
void AbsAlgorithm<L, D>::logfile_open() {
    logFile.open(full_result_path.c_str());
    if (logFile.fail()) {
        Print("!!! Cannot open experiment result file: ", full_result_path);
        exit(0);
    }
}

template<typename L, typename D>
void AbsAlgorithm<L, D>::result_log(){
    switch (stop_type) {
        case 1:
            logFile << std::setprecision(15) << iter << "\t" << grad_norm << "\t" << time << std::endl;
            break;
        case 2:
            logFile << std::setprecision(15) << iter << "\t" << fun_value << "\t" << time << std::endl;
            break;
        default:
            Print("There is no such stop type!");
            break;
    }
};

template<typename L, typename D>
void AbsAlgorithm<L, D>::Krcd_compute(L& K, const L& blocksize) {
    if(mu > 0)
        K = ceil(2 * exp(1) * data->n * (sqrt(1 + 1. / mu) - 1) / (D)blocksize + 1);
    else
        K = 10*data->n;
}

template<typename L, typename D>
void AbsAlgorithm<L, D>::primal_initial_set() {
    loss->w_initial(w, data, x, loss->b);
    grad_compute(grad,w, x);
    fun_value_compute(initial_fun_value,w, x);
    grad_norm_update(initial_grad_norm);
    grad_norm = initial_grad_norm;
    fun_value = initial_fun_value;
    switch (stop_type) {
        case 1:
            left = grad_norm;
            break;
        case 2:
            left = fun_value - fun_value_opt ;
            break;
        default:
            Print("There is no such stop type!");
            break;
    }
}



template<typename L, typename D>
void AbsAlgorithm<L, D>::check_blocksize(const L& blocksize) {
    if (blocksize > data->n) {
        Print("The size of block is larger than n! Please choose a new block size");
        exit(0);
    }
}

template<typename L, typename D>
void AbsAlgorithm<L, D>::reg_cord_cal(std::vector<D> &x,std::vector<D>& delta_x, const D& coefficient ,
                                         const std::vector<D> &cord_grad, const L& blocksize, const std::vector<L>& cord) {
    L tmp_cord;
    switch (reg_type) {
        case 1://indicator fun reg
            //Print("check 1");
            for (L row = 0; row < blocksize ; ++row) {
                tmp_cord = cord[row];
                delta_x[row] = -cord_grad[row] / (coefficient * v[tmp_cord]);
                if ((tmp_cord < data->nplus) && (x[tmp_cord] + delta_x[row] < 0)) {
                    delta_x[row] = -x[tmp_cord];
                    x[tmp_cord] = 0;
                } else {
                    x[tmp_cord] += delta_x[row];
                }
                //Print("cord", tmp_cord);
                //Print("delta", delta_x);
            }
            break;
        case 2://L1 reg
            D tmp_coef, tmp;
            //Print("check 2");
            for (L row = 0; row < blocksize ; ++row) {
                tmp_cord = cord[row];
                tmp_coef = coefficient * v[tmp_cord];
                tmp = tmp_coef * x[tmp_cord] - cord_grad[row];
                if (tmp < -1) {
                    delta_x[row] = (-cord_grad[row] + 1) / tmp_coef;
                    x[tmp_cord] = (tmp + 1) / tmp_coef;
                } else if (tmp > 1) {
                    delta_x[row] = (-cord_grad[row] - 1) / tmp_coef;
                    x[tmp_cord] = (tmp - 1) / tmp_coef;
                } else {
                    delta_x[row] = -x[tmp_cord];
                    x[tmp_cord] = 0;
                }
            }
            break;
        default:
            Print("There is no such reg type!");
            break;
    }
}

template<typename L, typename D>
void AbsAlgorithm<L, D>::grad_compute(std::vector<D>& grad, const std::vector<D>& w, const std::vector<D>& x){
    //Print("check,check");
    loss->cord_grad_update(grad,data,w,x,data->n,full_block);
    //Print("grad_before",grad);
    switch (reg_type) {
        case 1:
            for (L row = 0; row < data->nplus; ++row) {
                if ((x[row] == 0) && (grad[row] > 0))
                    grad[row] = 0;
            }
            break;
        case 2:
            for (L row = 0; row < data->n; ++row) {
                if (x[row] > 0)
                    grad[row] += reg_alpha;
                else if(x[row] < 0)
                    grad[row] -= reg_alpha;
                else{
                    if(grad[row] > 1)
                        grad[row] -= 1;
                    else if(grad[row] < -1)
                        grad[row] += 1;
                    else
                        grad[row] = 0;
                }
            }
            break;
        default:
            Print("There is no such reg type!");
            break;
    }
    // Print("grad",grad);
}

template<typename L, typename D>
void AbsAlgorithm<L, D>::grad_sum_compute(std::vector<D> &grad_sum, const std::vector<D> &w_y, const std::vector<D> &y,
                                          const std::vector<D> &w_z, const std::vector<D> &z, const D &theta_y) {

    L n = data->n;
    L m = data->m;
    std::vector<D> sum(n);
    std::vector<D> sum_w(m);
    for (L row = 0; row < n; ++row)
        sum[row] = theta_y * y[row] + z[row];
    for (L row = 0; row < m; ++row)
        sum_w[row] = theta_y * w_y[row] + w_z[row];
    grad_compute(grad_sum,sum_w,sum);

}

template<typename L, typename D>
void AbsAlgorithm<L, D>::fun_value_compute(D& fun_value,const std::vector<D>& w, const std::vector<D>& x){
    loss->fun_value_compute(fun_value,data,w,x);
    switch (reg_type) {
        case 1:
            break;
        case 2:
            for (L row = 0; row < data->n; ++row)
                fun_value += reg_alpha * abs(x[row]);
            break;
        default:
            Print("There is no such reg type!");
            break;
    }

};


#endif //ALGORITHM_ABSALGORITHM_H
