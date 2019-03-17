//
// Created by Lu Meng on 2018/11/27.
//
#include "LPuse.h"
int main() {
    std::string root_directory = "/Users/lumeng/Desktop/CCode/LP/data";
    const std::string all_experiments[10] = {"news20"};
    std::string experiment = all_experiments[0];
    std::vector<int> A_row_ptr;
    std::vector<int> A_col_idx;
    std::vector<double> A_value;
    std::vector<double> b;
    std::vector<double> c;

    /*
    int n,d;
    std::string full_data_path = root_directory + "/"+experiment+"/"+experiment;
    MatrixRead(full_data_path, n, d, A_row_ptr, A_col_idx, A_value,b);
    Print(class_sign<int,double>(b));
     */


    SVMtoLP<int,double>(root_directory,experiment);

    /*
    std::string full_path = root_directory +"/" + experiment +"/SVMtoLP";
    int mi,me,nb,nf;
    readMeta(full_path+"/meta",mi,me,nb,nf);
    int m = mi+me;
    int n = nb+nf;
    readMat(full_path +"/A", mi, n, A_row_ptr, A_col_idx, A_value);
    VectorRead(full_path +"/c",n,c);
    VectorRead(full_path +"/b",m,b);
    Print(m,n,A_row_ptr,A_col_idx,A_value);
     */




    /*
    int mi;
    int me;
    int nb;
    int nf;
    int n;
    std::vector<int> Ai_row_ptr;
    std::vector<int> Ai_col_idx;
    std::vector<double> Ai_value;
    root_directory = root_directory + "/" + experiment;
    readMeta(root_directory+"/meta", mi, me, nb, nf);
    Print("mi");
    Print(mi);
    Print("me");
    Print(me);
    Print("nb");
    Print(nb);
    Print("nf");
    Print(nf);
    n = nb+nf;
    readMat(root_directory+"/A", mi, n,Ai_row_ptr, Ai_col_idx, Ai_value);
    Print(Ai_row_ptr);
    Print(Ai_col_idx);
    Print(Ai_value);
    Print(mi,n,Ai_row_ptr,Ai_col_idx,Ai_value);
    */

    return 0;


}