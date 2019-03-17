// cdn_prox.cpp : Defines the entry point for the application.
//

#include "../use/useheader.h"
#include "L2LS.h"
#include "LS.h"



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

	std::string c_file_name = "c_4.txt";

	/*Print("Check L2LS_primal:");
	L2LS_primal<unsigned int,double> *inst = new L2LS_primal<unsigned int, double>(root_directory, experiment_name);
	Print("Check L2LS_dual:");
	L2LS_dual<unsigned int,double> *inst = new L2LS_dual<unsigned int, double>(root_directory, experiment_name);
	Print("Check LS_primal:");
	LS_primal<unsigned int, double> *inst = new LS_primal<unsigned int, double>(root_directory, experiment_name);
	*/
	Print("Check for LStoLP:");
	LPtoLS<unsigned int, double> *inst = new LPtoLS<unsigned int, double>(root_directory, experiment_name, c_file_name);


	unsigned int n = inst->n;
	unsigned int d = inst->d;
	Print("n:");
	Print(n);
	Print("d");
	Print(d);
	/*
	Print("lambda");
	Print(inst->lambda);
	 */
	Print("b:");
	Print(inst->b);
	Print("A");
	Print(n,d,inst->A_row_ptr,inst->A_col_idx,inst->A_value);
	Print("AT");
	Print(d,n,inst->AT_row_ptr,inst->AT_col_idx,inst->AT_value);
	Print("scalevalue:");
	Print(inst->scalevalue);
	Print("Diag:");
	Print(inst->ATL);

	std::vector<double> x(d,1);
	std::vector<double> w = Ax(n,inst->A_row_ptr,inst->A_col_idx,inst->A_value,x);
	w = VminusV(n,w,inst->b);
	Print("w:");
	Print(w);
	double primal;
	double dual;
	Print("PrimalDualGap:");
	Print(inst->PrimalDualGap(primal,dual,x,w));


	unsigned int tau  = d;
	std::vector<unsigned int> S(d);
	for (int i = 0; i < tau; ++i)
		S[i] = i;
	Print("Gs:");
	Print(inst->Gsub(tau,S,x,w));


	Print("Hs:");
	Print(tau,Msub(inst->AT_row_ptr, inst->AT_col_idx, inst->AT_value, inst->ATL, tau, S, inst->scalevalue));

	std::vector<double> delta_x(d,1);
	inst->w_update(w,tau,S,delta_x);
	Print("w");
	Print(w);
	delete inst;

//	std::cout << "\a\a\a\a\a\a\a" << std::endl;
	return 0;
}
