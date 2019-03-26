#include "../Problem/ProblemHeader.h"
#include "../LossFun/LossFunHeader.h"
#include "AbsAlgorithm.h"
#include "AbsAlgorithmSub.h"
#include "ParamPPA.h"
#include "PPA.h"
//#include "PCDM.h"
//#include "APPROX.h"

int main() {
    std::string root_directory = "/Users/lumeng/Desktop/CCode/LP/data";
    const std::string all_experiments[10] = {"Mat1000|1000_10000_0.0001_0.01", "Mat1000|1000_10000_0.0001_0.5", "test4", "test3",
                                             "Mat10000|10000_100000_0.0001_0.001", "Mat0|10000_60000_0.001_0.5",};
    std::string experiment = all_experiments[0];
    ProbData<unsigned int, long double> *inst = new ProbData<unsigned int, long double>(root_directory, experiment, "LPtoLP_inequ");
    //inst->data_print();
    //Print("b",inst->b);
    ParamPPA<long double>* param = new ParamPPA<long double>(1,0.9);

    /* AbsAlgorithm<unsigned int, long double>(data,sub_loss_type, sub_reg_type, sub_method_type,data_file_path_value)
     * sub_loss_type: 1. 0.5||Ae*x-be||^2 + 0.5||[Ai*x-bi]_{+}||^2 + 0.5*lambda||x-d||^2
     * sub_reg_type: 1. no regularization; 2. x[nb] >=0; 3. |x|
     * sub_method_type: 1. APPROX
     *
     * method->solver(x_initial, lambda_initial,epsilon_initial, delta_initial, sub_method_blocksize, max_iter, precision, result_iflog, method_type)
     * result_iflog: true; false
     * method_type: ""; "restart"
     */
    AbsAlgorithmSub<unsigned int, long double> *method = new PPA<unsigned int, long double>(inst, 1,1,1);
    param->varepsilon = 10;
    param->inv_idx_varepsilon = 1.1;
    param->inv_idx_delta = 1.1;
    param->exp_idx_theta = 2;
    param->theta = 1;
    unsigned int tau = 10;
    long double precision = 1e-6;
    std::vector<long double> x(inst->n, 0);
    std::vector<long double> lambda(inst->m, 0);
    method->solver(x,lambda,tau, 100, precision,param,false,"restart");

    method->result_print();

    delete method;
    delete param;
    delete inst;


    return 0;
}