#ifndef ALGORITHM_ABSALGORITHMSUB_H
#define ALGORITHM_ABSALGORITHMSUB_H

#include "AbsAlgorithm.h"
#include "APPROX.h"
#include "PCDM.h"

template<typename L, typename D>
class AbsAlgorithmSub:public AbsAlgorithm<L, D>{
protected:
    AbsAlgorithm<L, D>* sub_method;
    unsigned int sub_loss_type;
public:

    AbsAlgorithmSub(ProbData<L, D>* const data_inst, const unsigned int& sub_loss_type_value,
                    const unsigned int& sub_method_type);
    AbsAlgorithmSub(ProbData<L, D>* const data_inst, const unsigned int& sub_loss_type_value){}
    virtual ~AbsAlgorithmSub() {
        delete sub_method;
    }
};
template<typename L, typename D>
AbsAlgorithmSub<L, D>::AbsAlgorithmSub(ProbData<L, D>* const data_inst, const unsigned int& sub_loss_type_value,
                                       const unsigned int& sub_method_type):
        AbsAlgorithm<L, D>(data_inst,sub_loss_type), sub_loss_type(sub_loss_type_value) {
    switch (sub_method_type) {
        case 1:
            sub_method = new APPROX<L, D>(data_inst,sub_loss_type);
            break;
        case 2:
            sub_method = new PCDM<L, D>(data_inst,sub_loss_type);
            break;
        default:
            Print("There is no such sub_method type!");
            break;
    }
}


#endif //ALGORITHM_ABSALGORITHM_H
