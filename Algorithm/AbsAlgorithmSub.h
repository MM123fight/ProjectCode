#ifndef ALGORITHM_ABSALGORITHMSUB_H
#define ALGORITHM_ABSALGORITHMSUB_H

#include "AbsAlgorithm.h"
#include "APPROX.h"
#include "PCDM.h"
#include "ParamPPA.h"

template<typename L, typename D>
class AbsAlgorithmSub:public AbsAlgorithm<L, D>{
protected:
    AbsAlgorithm<L, D>* sub_method;
    unsigned int sub_loss_type;
public:
    AbsAlgorithmSub(ProbData<L, D>* const data_inst, const unsigned int& sub_loss_type,
    const unsigned int& sub_reg_type,const unsigned int& sub_method_type, const std::string& data_file_path="");
    virtual ~AbsAlgorithmSub() {
        delete sub_method;
    }
    virtual void solver(const std::vector<D>& x_initial,const std::vector<D>& lambda_initial,const L& sub_method_blocksize,
                        const L& max_iter, const D& precision, ParamPPA<D>* const param,
                        const bool& result_iflog = false, const std::string& method_type = ""){}
};
template<typename L, typename D>
AbsAlgorithmSub<L, D>::AbsAlgorithmSub(ProbData<L, D>* const data_inst, const unsigned int& sub_loss_type,
                                       const unsigned int& sub_reg_type,const unsigned int& sub_method_type,
                                       const std::string& data_file_path):
        AbsAlgorithm<L, D>(data_inst,sub_loss_type,sub_reg_type,data_file_path),sub_loss_type(sub_loss_type){
    switch (sub_method_type) {
        case 1:
            sub_method = new APPROX<L, D>(data_inst,sub_loss_type,sub_reg_type);
            break;
        default:
            Print("There is no such sub_method type!");
            break;
    }
}


#endif //ALGORITHM_ABSALGORITHM_H
