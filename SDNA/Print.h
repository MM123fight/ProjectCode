//
// Created by Lu Meng on 2018/5/24.
//

#ifndef SDNA_PRINT_H
#define SDNA_PRINT_H

#include "../use/stdafx.h"
template <typename D>
void ReportError(std::ofstream &logFile, const std::vector<D> &b){
    if (b.size() == 0) {
        std::cout << "!!! The data is not loaded! Please check back gain." << std::endl;
        exit(0);
    }
    if (!logFile.is_open()) {
        std::cout << "!!!The result file is not open. The data will not be recorded." << std::endl;
        exit(0);
    }
}

template < typename L,typename D>
void PrintInital(const L &tau, const D &fun_val, const L &precision, std::ofstream &logFile) {
    std::cout << "The initial function value:" << fun_val << ";tau =" << tau << std::endl;
    logFile << std::setprecision(precision) << 1e-15 << "\t" << log10(fun_val) << "\t" << 1e-15 << std::endl;
}

template <typename L, typename D>
void PrintData(const L &itertau, const int &precision, const D &fun_val,
                   std::ofstream &logFile, const D &elapsedTime, const bool &screenprint) {
    if(screenprint == true) {
        std::cout << itertau << "\tFunction Value: " << fun_val << ";\tElapsed time: " <<
                  elapsedTime << " s" << std::endl;
    }
    logFile << std::setprecision(precision) << itertau << "\t" << fun_val << "\t" << elapsedTime << std::endl;
}

template < typename L,typename D>
void PrintInital(const L &tau, const D &gap, const D &primal, const D &dual,
                 const L &precision, std::ofstream &logFile) {
    std::cout << "0\tDuality Gap: " << gap << " = " << primal << " - " << dual << "; tau = " << tau << std::endl;
    logFile << std::setprecision(precision) << 1e-15 << "\t" << gap << "\t" <<
            primal << "\t" << dual << "\t" << 1e-15 <<"\t"<<1e-15<< std::endl;
}

template <typename L, typename D>
void PrintData(const L &itertau, const int &precision, const D &gap, const D &primal, const D&dual,
               std::ofstream &logFile, const D &elapsedTime, const D& elapsedTimeSub, const bool &screenprint) {
    if(screenprint == true) {
        std::cout << itertau << "\t\tDuality Gap: " << gap << " = " << primal
                  << " - " << dual << ";\tElapsed time: " << elapsedTime << " s" << std::endl;
    }
    logFile << std::setprecision(precision) << itertau << "\t" <<gap << "\t" <<
            primal << "\t" << dual << "\t" << elapsedTime << "\t" << elapsedTimeSub << std::endl;
}


template < typename L, typename D>
void PrintFinal(const D &elapsedTime, const D &elapsedTimeSampling, const D &elapsedTimeComputingHessian,
                const D &elapsedTimeCD, const D &elapsedTimeLbfgsb, const D &skipped_cd, const L &skipped_lbfgsb,
                const L &itertau, std::ofstream &logFile) {
    std::cout << "Finished." << std::endl;
    std::cout << "Total epochs: " << itertau << std::endl;
    std::cout << "Elapsed time: " << elapsedTime << " = " <<
              "Sampling (" << elapsedTimeSampling << ") + Hessian(" << elapsedTimeComputingHessian <<
              ") +CD (" << elapsedTimeCD << ") + L-BFGS-B (" <<
              elapsedTimeLbfgsb << ")" << std::endl;
    std::cout << "Percentage of CD iterations skipped: " << 100 * skipped_cd / (D) itertau << " %" << std::endl;
    std::cout << "Percentage of LBFGS iterations skipped: " << 100 * skipped_lbfgsb / (D) itertau << " %" << std::endl;

    logFile << "#Elapsed time: " << elapsedTime << " = " <<
                           "Sampling (" << elapsedTimeSampling << ") + Hessian(" << elapsedTimeComputingHessian <<
                           ") +CD (" << elapsedTimeCD << ") + L-BFGS-B (" <<
                           elapsedTimeLbfgsb << ")" << std::endl;

}
template < typename L, typename D>
void PrintFinal(const D &elapsedTime, const D &elapsedTimeSampling, const D &elapsedTimeSub,
                const L &itertau) {
    std::cout << "Finished." << std::endl;
    std::cout << "Elapsed time: " << elapsedTime << " = " <<
              "Sampling (" << elapsedTimeSampling << ") + sub (" <<
              elapsedTimeSub << ")" << std::endl;
}

template < typename L, typename D>
void PrintFinal(const D &elapsedTime, const D &elapsedTimeSampling, const D &elapsedTimeComputingHessian, const D &elapsedTimeSub,
                const L &itertau) {
    std::cout << "Finished." << std::endl;
    std::cout << "Elapsed time: " << elapsedTime << " = " <<
              "Sampling (" << elapsedTimeSampling << ")+Hessian(" << elapsedTimeComputingHessian
              << ") + sub (" << elapsedTimeSub << ")" << std::endl;
}

#endif //SDNA_PRINT_H
