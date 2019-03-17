// cdn_prox.cpp : Defines the entry point for the application.
//

#include "../use/stdafx.h"
#include "../use/Print.h"
#include "LeastSquare.h"
#include "HingeLoss.h"
#include "L2LS.h"



int main() {
	/* This makes cin and cout much faster (more than 10x) at the cost of
	 * disabling usage of IO functions in c, e.g. scanf/printf.
	 * Check the details on the following web:
	 * https://blog.csdn.net/qq_33248299/article/details/52144485
	 */
	std::ios_base::sync_with_stdio(false);

	// The data and code directory;
    std::string root_directory = "/Users/lumeng/Desktop/CCode/data";

    //In order to obtain the random number based on time
	srand((unsigned)time(NULL));

    //data file names
	const std::string all_experiments[10] = {"LS5_4.txt","gisette_scale","madelon.txt","w2a"};
	// modify here to change experiment and you can also input the file name by yourself
	std::string experiment_name = all_experiments[0];

	// Check the GenInfo.h
	GenInfo<int,double> *inst = new L2Loss_primal<int, double>(root_directory, experiment_name);
	int n = inst->n;
	int d = inst->d;
	Print(n);
	Print(d);
	Print(inst->b);
	Print("A");
	Print(n,d,inst->A_row_ptr,inst->A_col_idx,inst->A_value);
	Print("AT");
	Print(d,n,inst->AT_row_ptr,inst->AT_col_idx,inst->AT_value);


	/*

	std::vector<std::vector<double>> A = std::vector<std::vector<double>>(n);
	for (int i = 0; i < n; ++i) {
		A[i].resize(m);
		for (int j = inst->A_csr_row_ptr[i]; j < inst->A_csr_row_ptr[i+1]; ++j) {
			A[i][inst->A_csr_col_idx[j]]= inst->A_csr_values[j];
		}
		Print(A[i]);
	}
	Print();

	std::vector<std::vector<double>> AT = std::vector<std::vector<double>>(m);
	for (int i = 0; i < m; ++i) {
		AT[i].resize(n);
		for (int j = inst->A_csc_col_ptr[i]; j < inst->A_csc_col_ptr[i+1]; ++j) {
			AT[i][inst->A_csc_row_idx[j]] = inst->A_csc_values[j];
		}
		Print(AT[i]);
	}
	Print();

	std::vector<std::vector<double>> ATA = std::vector<std::vector<double>>(m);
	for (int i = 0; i < m; ++i) {
		ATA[i].resize(m);
		for (int j = 0; j < m; ++j) {
			ATA[i][j] = inst->vector_multiple_primal(i,j);
		}
		Print(ATA[i]);
	}
	Print();

	std::vector<std::vector<double>> AAT = std::vector<std::vector<double>>(n);
	for (int i = 0; i < n; ++i) {
		AAT[i].resize(n);
		for (int j = 0; j < n; ++j) {
			AAT[i][j] = inst->vector_multiple_dual(i,j);
		}
		Print(AAT[i]);
	}
	Print();


	//! I find that if every time we run the algorithm from the beginning,
	//! we will generate the same samples.(Is it good?)
	int tau = 3;
	std::vector<int> S;
	std::vector<double> Qdata;
	param.solve_from_dual = true;
	gsl_rng *gsl_rng_r = gsl_rng_alloc(gsl_rng_default);
	for (int i = 0; i < 10; ++i) {
		inst->BlockGen(tau,S,param.solve_from_dual,gsl_rng_r);
		std::cout << "The chosen sample:";
		Print(S);
		std::cout << "The corresponding Hessian:\n";
		inst->BlockHessian(tau, Qdata, S, param.solve_from_dual);
		for (int i = 0; i < tau; ++i) {
			for (int j = 0; j < tau; ++j) {
				std::cout << Qdata[i * tau + j] << " ";
			}
			Print();
		}
		Print();
	}
	gsl_rng_free(gsl_rng_r);
	*/
	/*
	int n =100;
	std::cout << 1./100 << std::endl;
	L2Loss<int, double> *inst = new L2Loss<int,double>(root_directory, experiment_name);
	std::string row = to_string(inst->n);
	std::string col = to_string(inst->m);
	std::string result_name_begin = "result3/LS" + row + "_" + col + "_";
	std::string result_name;
	std::vector<int> tau;
	tau.push_back(2);
	Print(inst->Li);
	//Find why inst->setLambda(10.0) will make the code wrong?
	//inst->setLambda(10.0);
	L2param.per_output = 2*inst->n;
	for (int i = 0; i < tau.size(); ++i) {
		std::string tau_name = to_string(tau[i]);
		result_name = result_name_begin + tau_name + "_SDNA.txt";
		inst->setResultFile(result_name);
		inst->SparseSDNA(tau[i]);
	}
	 */


////	HingeLoss<int, double> *inst = new HingeLoss<int,double>(root_directory, experiment_name);
////
////
////
////	int tau = 5;
////
////	std::string row = to_string(inst->m);
////	std::string col = to_string(inst->n);
////	std::string tau_name = to_string(tau);
////	std::string result_name;
////
////	result_name = "hingeloss"+  row + "_" + col + "_" + tau_name + ".txt";
////
////	inst->setResultFile(result_name);
////	inst->SparseSDNA(tau);
//
//
////	LeastSquare<int, double> *inst = new LeastSquare<int,double>(root_directory, experiment_name, "c_4.txt");
//
//
//	 std::vector<int> tau;
//    for (int i = 2; i < 30; ++i) {
//        tau.push_back(i);
//    }
//
////     std::cout << "n:" << std::endl;
////     Print(inst->n);
////     std::cout << "Li:" << std::endl;
////     Print(inst->Li);
//
//
//
//
//
////     inst->setBound_idx(0.5 * inst->m);
//     std::string row = to_string(inst->n);
//     std::string col = to_string(inst->m);
//     std::string result_name_begin = "result3/LS" + row + "_" + col + "_";
//     std::string result_name;
//
//     L2param.max_time = 1e4;
//	 L2param.tol = 1e-10;
//	 L2param.per_output =1000;
//	 L2param.screenprint= false;
//	 L2param.theta = 1e-9;
//
//	std::string tau_name;
////	result_name = result_name_begin +  "_1.txt";
////	inst->setResultFile(result_name);
////	inst->SparsePCDM(1);
//	for (int i = 0; i < tau.size(); ++i) {
//
//        tau_name = to_string(tau[i]);
////		result_name = result_name_begin + tau_name + "_PCDM.txt";
////		inst->setResultFile(result_name);
////		inst->SparsePCDM(tau[i]);
//
//        result_name = result_name_begin + tau_name +"_SDNA.txt";
//        inst->setResultFile(result_name);
//        inst->SparseSDNA(tau[i]);
//	}
	delete inst;
//	std::cout << "\a\a\a\a\a\a\a" << std::endl;
	return 0;
}