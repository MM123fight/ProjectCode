// cdn_prox.cpp : Defines the entry point for the application.
//

#include "../use/useheader.h"
#include "GenInfo.h"



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
	const std::string all_experiments[10] = {"LS1000_10_0.3_0.5","LS10000_1000_0.3_0.0001","LS100000_1000_0.3_0.0001"};
	std::string experiment = all_experiments[0];

	root_directory = root_directory+ "/" + experiment;
	std::string experiment_name =  experiment + ".txt";
	//experiment_name = "scaled_" + experiment_name;
	std::string file = "scaled.txt";
	GenInfo<unsigned int, double> *inst = new GenInfo<unsigned int, double>(root_directory, experiment_name, false);
	std::vector<double> Leverage;
	std::string Leverage_path = root_directory + "/leverage_" + experiment_name;
	VectorRead(Leverage_path, inst->n, Leverage);
	double epsilon = (1 - param.leverage_approx) / (1 + param.leverage_approx);
	Print("alpha");
	Print(1/sqrt(epsilon));
	double c = -3 * log(1 - param.leverage_prob) / log(inst->d);
	Print("c");
	Print(c);
	Print("log d");
	Print(log(inst->d));
	Print("coefficient");
	Print(c/sqrt(epsilon)*log(inst->d));
	Print("possible leverage:");
	Print(sqrt(epsilon)/(log(inst->d)*c));
	int tau;
	std::vector<unsigned int> S;
	for ( int i = 0; i < 100; ++i) {
		leverage_sample(S, inst->n, inst->d, Leverage, epsilon, c);
		tau = S.size();
		Print(S);
	}
	//std::vector<double> x = inst->ExpectMatrixInverseEinv(root_directory, experiment_name);
	//Print(inst->n,x);

	delete inst;
//	std::cout << "\a\a\a\a\a\a\a" << std::endl;
	return 0;
}
