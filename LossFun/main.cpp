//
// Created by Lu Meng on 2018/11/27.
//
#include "../Problem/ProblemHeader.h"
#include "plus_squared_loss.h"

int main() {
    std::string root_directory = "/Users/lumeng/Desktop/CCode/LP/data";
    const std::string all_experiments[10] = {"LS6_4","test3","Mat0|1000_10000_0.01_0.5","Mat0|10000_60000_0.001_0.5"};
    //std::string experiment = all_experiments[0] + "/SVMtoLP";
    std::string experiment = all_experiments[1];

    Print("test");
    ProbData<unsigned int, double> *inst = new ProbData<unsigned int, double>(root_directory, experiment,"LP");

    double alpha = 1;
    unsigned int m = inst->m;
    unsigned int n = inst->n;

    AbsLoss<unsigned int, double>* loss;
    loss = new plus_squared_loss<unsigned int, double>(alpha,inst->b,inst->c);
    Print("test");

    inst->data_print();

    std::vector<double> w(m);
    std::vector<double> x(n,3);
    Print("x",x);
    std::vector<double> y(n,1);
    std::vector<double> z(n,2);
    std::vector<double> w_z(m);
    std::vector<double> w_y(m);
    loss->w_initial(w,inst,x,inst->b);
    Print("w",w);
    loss->w_initial(w_z,inst,z,inst->b);
    Print("w_z",w_z);
    loss->w_initial(w_y,inst,y);
    Print("w_y",w_y);
    std::vector<double> grad(n);
    loss->grad_compute(grad, inst, w, x);
    Print("grad",grad);
    loss->grad_sum_compute(grad,inst,w_y,y,w_z,z,1);
    Print("grad",grad);


    delete inst;
    delete loss;
    return 0;
}