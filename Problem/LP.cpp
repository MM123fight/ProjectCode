//
// Created by Lu Meng on 2018/11/27.
//


#include "ProbData.h"
int main() {
    std::string root_directory = "/Users/lumeng/Desktop/CCode/LP/data";
    const std::string all_experiments[10] = {"test3","Mat100|100_1000_0.01_0.5","Mat0|1000_10000_0.01_0.5","Mat0|10000_60000_0.001_0.5"};
    //std::string experiment = all_experiments[0] + "/SVMtoLP";
    std::string experiment = all_experiments[1];
    const std::string data_type[10] = {"LPtoLP_inequ","LP", "LPtoLP", "LPtoLS", "LPtoLS_equ","LS", "LS_equ"};
    std::string data_tp = data_type[1];

    ProbData<unsigned int, double> *inst = new ProbData<unsigned  int, double>(root_directory, experiment, data_tp);
    Print(inst->data_type);
    //inst->data_print();
    Print("b.size", inst->b.size());
    Print("b", inst->b);
    unsigned int nb = inst->nb;
    unsigned int m = inst->m;

    data_tp = data_type[0];
    inst = new ProbData<unsigned  int, double>(root_directory, experiment, data_tp);
    Print(inst->data_type);
    //inst->data_print();
    Print("b.size", inst->b.size());
    for (int i = 0; i < inst->m; ++i)
        std::cout << inst->b[i + nb] << "\t";

    /*
    data_tp = data_type[2];
    inst = new ProbData<unsigned  int, double>(root_directory, experiment, data_tp);
    Print(inst->data_type);
    inst->set_mplus(1);
    inst->set_nplus(2);
    Print(inst->mplus);
    Print(inst->nplus);
    inst->data_print();
     */

    delete inst;
    return 0;
}