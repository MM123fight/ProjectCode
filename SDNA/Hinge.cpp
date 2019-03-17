// cdn_prox.cpp : Defines the entry point for the application.
//

#include "../use/useheader.h"
#include "Hinge.h"



int main() {
	/* This makes cin and cout much faster (more than 10x) at the cost of
	 * disabling usage of IO functions in c, e.g. scanf/printf.
	 * Check the details on the following web:
	 * https://blog.csdn.net/qq_33248299/article/details/52144485
	 */
	std::ios_base::sync_with_stdio(false);
	// The data and code directory;
    std::string root_directory = "/Users/lumeng/Desktop/CCode/data/Hinge";
    //In order to obtain the random number based on time
	srand((unsigned)time(NULL));

    //data file names
	const std::string all_experiments[10] = {"LS5_4","LS100000_100_0.3_0.0001","LS10000_1000_0.3_0.0001","LS100000_1000_0.3_0.0001"};
	std::string experiment = all_experiments[1];

	root_directory = root_directory+ "/" + experiment;
	std::string experiment_name =  experiment + ".txt";
	//experiment_name = "scaled_" + experiment_name;
	std::string file = "scaled.txt";

	param.screenprint = false;
	param.per_output = 0.1;
	param.max_time = 200;
	std::vector<double> multiple;
	//multiple= {0.0001,0.001,0.01,0.1,1,10,100,1000};
    multiple = {1};
	std::string result_name;
	std::vector<double> x;

	std::string method_types[10] = {"SDCA","newSDCA","d_APPROX","dSDNA"};
	std::vector<unsigned int> tau_choices;
	tau_choices= {1,2,4,8,16,32,64};
    //tau_choices = {1};
	std::string method = method_types[2];
	unsigned int tau;

	if(method == "dSDNA") {
		Print("Check Hinge_dSDNA:");
        Hinge_dual<unsigned int, double> *inst = new Hinge_dual<unsigned int, double>(root_directory, experiment_name, false);
        //Print(inst->n,inst->d,inst->A_row_ptr,inst->A_col_idx,inst->A_value);
        //Print(inst->lower_bound);
		//Print(inst->upper_bound);
		std::string result_dSDNA;
		for (int j = 0; j < tau_choices.size(); ++j) {
			tau = tau_choices[j];
			result_dSDNA = "dSDNA_size" + to_string(tau) + "_lambda";
			for (int i = 0; i < multiple.size(); ++i) {
				param.lambda_multiple = multiple[i];
				result_name = result_dSDNA + to_string(param.lambda_multiple) + "_" + file;
				inst->setResultFile(result_name);
				x = inst->SDNA(tau, "GM");
				Print(x);
			}
		}
		delete inst;
		inst = NULL;
	}else if(method == "newSDCA"){
		Print("Check Hinge_newSDCA");
        Hinge_dual<unsigned int, double> *inst = new Hinge_dual<unsigned int, double>(root_directory,
																						experiment_name, false);

		std::string result_SDCAnew;
		for (int j = 0; j < tau_choices.size(); ++j) {
			tau = tau_choices[j];
			result_SDCAnew = "SDCAnew_size" + to_string(tau) + "_lambda";
			for (int i = 0; i < multiple.size(); ++i) {
				param.lambda_multiple = multiple[i];
				result_name = result_SDCAnew + to_string(param.lambda_multiple) + "_" + file;
				inst->setResultFile(result_name);
				x = inst->PCDM(tau, "newPCDM");
			}
		}
		delete inst;
		inst = NULL;
	}else if(method == "SDCA"){
		Print("Check Hinge_SDCA");
        Hinge_dual<unsigned int, double> *inst = new Hinge_dual<unsigned int, double>(root_directory,
																						experiment_name, false);
		std::string result_SDCA;
		for (int j = 0; j < tau_choices.size(); ++j) {
			tau = tau_choices[j];
			result_SDCA= "SDCA_size" + to_string(tau) + "_lambda";
			for (int i = 0; i < multiple.size(); ++i) {
				param.lambda_multiple = multiple[i];
				result_name = result_SDCA + to_string(param.lambda_multiple) + "_" + file;
				inst->setResultFile(result_name);
				x = inst->PCDM(tau, "PCDM");
			}
		}
		delete inst;
		inst  = NULL;

	}else if(method == "d_APPROX"){
		Print("Check Hinge_d_Approx:");
        Hinge_dual<unsigned int, double> *inst = new Hinge_dual<unsigned int, double>(root_directory,
																						experiment_name, false);
		std::string result_d_Approx;
		for (int j = 0; j < tau_choices.size(); ++j) {
			tau = tau_choices[j];
			result_d_Approx = "d_Approx_size" + to_string(tau) + "_lambda";
			for (int i = 0; i < multiple.size(); ++i) {
				param.lambda_multiple = multiple[i];
				result_name = result_d_Approx + to_string(param.lambda_multiple) + "_" + file;
				inst->setResultFile(result_name);
				x = inst->APPROX(tau);
			}
		}
		delete inst;
		inst = NULL;
	}
	else{
		Print("No such method");
	}


//	std::cout << "\a\a\a\a\a\a\a" << std::endl;
	return 0;
}
