#include "../Problem/ProblemHeader.h"
#include "../LossFun/LossFunHeader.h"
#include "AbsAlgorithm.h"
#include "ParamPPA.h"
//#include "PPA.h"
#include "PCDM.h"
#include "APPROX.h"

int main() {
    /* AbsAlgorithm<unsigned int, long double>(data,loss_type, reg_type, data_file_path_value)
     * loss_type: 1. 0.5||Ae*x-be||^2 + 0.5||[Ai*x-bi]_{+}||^2 + 0.5*lambda||x-d||^2
     * reg_type: 1. no regularization; 2. x[nb] >=0; 3. |x|
     *
     * method->solver(x_initial, blocksize, max_iter, precision, stop_type = 1, result_iflog, method_type)
     * stop_type: 1. grad_norm; 2. fun_value - fun_opt_value
     * result_iflog: true; false
     * method_type: ""; "restart"
     */
    std::string root_directory = "/Users/lumeng/Desktop/CCode/LP/data";
    const std::string all_experiments[10] = {"Mat1000|1000_10000_0.0001_0.01", "Mat10000|10000_100000_0.0001_0.001", "test1", "test3",
                                             "Mat0|1000_10000_0.01_0.5", "Mat0|10000_60000_0.001_0.5"};
    //std::string experiment = all_experiments[0] + "/SVMtoLP";
    std::string experiment = all_experiments[0];
    //std::string full_result_path = root_directory + "/" + experiment + "/result/test.txt";
    ProbData<unsigned int, long double> *inst = new ProbData<unsigned int, long double>(root_directory, experiment, "LPtoLP_inequ");
    //Print("inst->b", inst->b);
    Print("mplus", inst->mplus);
    Print("nplus", inst->nplus);

    std::vector<long double> initial_x;
    std::vector<long double> sub_b;
    std::vector<long double> sub_d;
    std::string sub_b_path = inst->data_file_path + "/result/sub_b";
    VectorRead(sub_b_path,inst->m,sub_b);
    std::string sub_d_path = inst->data_file_path + "/result/sub_d";
    VectorRead(sub_d_path,inst->n,sub_d);
    std::string initial_x_path = inst->data_file_path + "/result/initial_x";
    VectorRead(initial_x_path,inst->n,initial_x);

    AbsAlgorithm<unsigned int, long double> *method = new APPROX<unsigned int, long double>(inst, 1, 1);
    long double sigma = 4;
    long double alpha = 1./(sigma*sigma);
    unsigned int tau = 1;
    long double precision = 1e-15;
    method->set_loss(alpha,sub_b,sub_d);
    method->solver(initial_x, tau, 100, precision, 1,true,"restart");
    //method->restart_step();
    method->result_print();
    /*
    x.clear();
    x.resize(inst->n,0);
    tau = 1;
    method->solver(x,tau , 10000, precision, 1, false,"restart");

    method->initial_print();
    method->result_print();
     */
    /*
    method->step();
    method->initial_print();
    method->result_print();
     */
    //inst->data_print();

    Print("\n");
    Print("x", method->x);

    /*
    double value = 0;
    unsigned int idx;
    Print("size", inst->m);
    Print("nb",inst->nb);
    Print("nf",inst->nf);
    Print("size", inst->A_row_ptr[inst->m]);
    Print("size", inst->A_row_ptr[inst->m-1]);
    Print("A_row_ptr", inst->A_row_ptr);
    for (int col = inst->A_row_ptr[inst->m-1]; col < inst->A_row_ptr[inst->m]; ++col) {
        idx = inst->A_col_idx[col];
        if((idx < inst->nb)||(idx >= inst->n-inst->nf))
            value += inst->A_value[col] * method->x[idx];
    }
    Print("LP objective value", value);
     */

    /*
    method->step();

    method->initial_print();
    method->result_print();
    Print("\n");
    method->solver(method->x,tau,2,precision,2,true);
    method->initial_print();
    method->result_print();
     */






    delete method;
    delete inst;



    return 0;
}