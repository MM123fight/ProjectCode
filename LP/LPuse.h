//
// Created by Lu Meng on 2018/11/27.
//

#ifndef SDNA_LPUSE_H
#define SDNA_LPUSE_H

#include "../use/useheader.h"

template <typename L>
void readMeta(const std::string &full_data_path ,L &mi, L &me, L &nb, L &nf);

template <typename L>
void writeMeta(const std::string &full_data_path ,const L &mi, const L &me, const L &nb, const L &nf);

template <typename L, typename D>
void readMat(const std::string &full_data_path, const L &m, const L &n, std::vector<L> &A_row_ptr,
             std::vector<L> &A_col_idx, std::vector<D> &A_value );
template <typename L, typename D>
void writeMat(const std::string &full_data_path, const L &m, const L &n, const std::vector<L> &A_row_ptr,
              const std::vector<L> &A_col_idx, const std::vector<D> &A_value );
template <typename L, typename D>
L check_class(const D &y, const L &k, const std::vector<D> &sign);

template <typename L, typename D>
std::vector<D> class_sign( const std::vector<D> &y);

template <typename L, typename D>
void SVMtoLP(const std::string &root_directory, const std::string &experiment, D lambda = 1);



template <typename L>
void readMeta(const std::string &full_data_path ,L &mi, L &me, L &nb, L &nf){
    std::ifstream fin(full_data_path.c_str());
    if (fin.fail()) {
        std::cerr << "Cannot open data file: " << full_data_path << std::endl;
        exit(0);
    }
    std::string tmp;
    fin >> tmp >> nb;
    fin >> tmp >> nf;
    fin >> tmp >> mi;
    fin >> tmp >> me;
    fin.close();
}

template <typename L>
void writeMeta(const std::string &full_data_path ,const L &mi, const L &me, const L &nb, const L &nf){
    std::ofstream fout;
    fout.open(full_data_path.c_str());
    if (fout.fail()) {
        std::cerr << "Cannot open data file: " << full_data_path << std::endl;
        exit(0);
    }

    fout << "nb" << " "<< nb << "\n";
    fout << "nf" << " "<< nf << "\n";
    fout << "mi" << " "<< mi << "\n";
    fout << "me" << " "<< me;
    fout.close();
}

template <typename L, typename D>
void readMat(const std::string &full_data_path, const L &m, const L &n, std::vector<L> &A_row_ptr,
             std::vector<L> &A_col_idx, std::vector<D> &A_value ){

    L i,j;
    D val;
    std::ifstream fin(full_data_path.c_str());
    if( fin.fail() ){
        std::cerr << "fail to open " << full_data_path << std::endl;
        exit(0);
    }
    fin >> i >> j >> val; //filter one line
    if( i != m || j != n ){
        std::cerr << "dimension in " << full_data_path << " does not match that in meta file" << std::endl;
        exit(0);
    }
    L row_idx = 1;
    L tmp = 0;
    A_row_ptr.push_back(0);
    while( !fin.eof() ){
        fin >> i >> j >> val;

        if( i-1 >= m || j-1 >= n ){
            std::cerr << "index:" << "(" << i-1 << ", " << j-1 << ") out of bound when reading " << full_data_path << std::endl;
            exit(0);
        }
        A_col_idx.push_back(j-1);
        A_value.push_back(val);
        if(i != row_idx){
            A_row_ptr.push_back(tmp);
            ++row_idx;
        }
        ++tmp;
    }
    if(tmp != 0)
        A_row_ptr.push_back(tmp);
    fin.close();
}

template <typename L, typename D>
void writeMat(const std::string &full_data_path, const L &m, const L &n, const std::vector<L> &A_row_ptr,
             const std::vector<L> &A_col_idx, const std::vector<D> &A_value ){

    L i,j;
    D val;
    std::ofstream fout;
    fout.open(full_data_path.c_str());
    if( fout.fail() ){
        std::cerr << "fail to open " << full_data_path << std::endl;
        exit(0);
    }
    fout << m << " " << n << " " << "0.0" << "\n"; //filter one line
    for (L row = 0; row < m; ++row) {
        for (L col = A_row_ptr[row]; col < A_row_ptr[row+1]; ++col)
            fout << row+1 << " " << A_col_idx[col] + 1 << " " << A_value[col] << "\n";
    }
    fout.close();
}

template <typename L>
void writeMat(const std::string &full_data_path, const L &n){

    std::ofstream fout;
    fout.open(full_data_path.c_str());
    if( fout.fail() ){
        std::cerr << "fail to open " << full_data_path << std::endl;
        exit(0);
    }
    fout  << "0 " << n << " " << "0.0" << "\n"; //filter one line

    fout.close();
}
template <typename L, typename D>
void Libsvm_writeMat(const std::string &full_data_path, const L &m, const L &n, const std::vector<L> &A_row_ptr,
                     const std::vector<L> &A_col_idx, const std::vector<D> &A_value, const std::vector<D> &b ) {
    std::ofstream fout;
    fout.open(full_data_path.c_str());
    for (L i = 0; i < m; i++) {
        fout << b[i] << " ";
        for (L j = A_row_ptr[i]; j < A_row_ptr[i+1]; j++)
            fout << A_col_idx[j] + 1 << ":" << A_value[j] << " ";
        fout << std::endl;
    }
    fout.close();
}




template <typename L, typename D>
L check_class(const D &y, const L &k, const std::vector<D> &sign){
    if(sign.size()!=k){
        Print("The size of class is not accordant");
        exit(0);
    }
    L begin = 0;
    L end = k -1;
    L check_point = (begin+end)/2;
    L class_type;
    while(check_point!=begin){
        if(y == sign[check_point]) {
            class_type = check_point + 1;
            break;
        }
        else if(y < sign[check_point])
            end = check_point;
        else
            begin = check_point;
        check_point = (begin+end)/2;
    }
    if(y==sign[check_point]){
        class_type = check_point+1;
    }else if (y==sign[end]){
        class_type = end+1;
    }else{
        Print("It does not belong to any class");
        exit(0);
    }

    return class_type;
};

template <typename L, typename D>
std::vector<D> class_sign( const std::vector<D> &y){
    std::vector<D> sign;
    L iter = 0;
    sign.push_back(y[0]);
    L sign_size = 1;
    for (L i = 1; i < y.size(); ++i) {
        if(y[i] > sign[iter]){
            ++iter;
            while((y[i] > sign[iter])&&(iter<sign_size))
                ++iter;
            if(((iter<sign_size)&&(y[i]!= sign[iter]))||(iter == sign_size)) {
                sign.insert(sign.begin() + iter, y[i]);
                ++sign_size;
            }
        }else if(y[i] < sign[iter]){
            --iter;
            while((y[i] < sign[iter])&&(iter>=0))
                --iter;
            if(((iter>=0)&&(y[i]!= sign[iter]))||(iter == -1)) {
                sign.insert(sign.begin() + iter + 1, y[i]);
                ++sign_size;
            }
            if(iter == -1)
                iter = 0;
        }
    }
    return sign;
};

template <typename L, typename D>
void SVMtoLP(const std::string &root_directory, const std::string &experiment, D lambda){
    std::string full_data_path = root_directory + "/" + experiment + "/" + experiment;
    std::vector<L> X_row_ptr;
    std::vector<L> X_col_idx;
    std::vector<D> X_value;
    std::vector<D> ne_X_value;
    std::vector<D> y;
    std::vector<D> sign;
    L N,d;
    MatrixRead(full_data_path, N, d, X_row_ptr, X_col_idx, X_value, y);
    ne_X_value.resize(X_value.size());
    for (L row = 0; row < X_value.size(); ++row)
        ne_X_value[row] = - X_value[row];
    sign = class_sign<L,D>(y);
    L k = sign.size();
    L m = (k-1)*N;
    L n = 2*k*d+N;
    std::vector<D> c(n,1);
    for (L row = 0; row < 2*k*d; ++row)
        c[row] = lambda;
    std::vector<D> b(m,-1);
    std::vector<L> A_row_ptr( m+1,0);
    std::vector<L> A_col_idx;
    std::vector<D> A_value;


    L begini, endi, beginj,endj;
    L size, classyi;

    for (L i = 0; i < N; ++i) {
        classyi = check_class(y[i], k, sign);
        begini = X_row_ptr[i];
        endi = X_row_ptr[i+1];
        size = 4*(endi - begini)+1;
        for (L j = 0; j < classyi-1; ++j){
            A_row_ptr[i*(k-1)+j+1] = A_row_ptr[i*(k-1)+j]+size;
            A_value.insert(A_value.end(),X_value.begin()+begini, X_value.begin()+endi);
            A_value.insert(A_value.end(),ne_X_value.begin()+begini, ne_X_value.begin()+endi);
            A_value.insert(A_value.end(),ne_X_value.begin()+begini, ne_X_value.begin()+endi);
            A_value.insert(A_value.end(),X_value.begin()+begini, X_value.begin()+endi);
            A_value.push_back(-1);
            for (L row = begini; row < endi; ++row)
                A_col_idx.push_back(X_col_idx[row]+2*j*d);
            for (L row = begini; row < endi; ++row)
                A_col_idx.push_back(X_col_idx[row]+(2*j+1)*d);
            for (L row = begini; row < endi; ++row)
                A_col_idx.push_back(X_col_idx[row]+2*(classyi-1)*d);
            for (L row = begini; row < endi; ++row)
                A_col_idx.push_back(X_col_idx[row]+(2*classyi-1)*d);
            A_col_idx.push_back(2*k*d+i);
        }
        for (L j = classyi; j < k; ++j){
            A_row_ptr[i*(k-1)+j] = A_row_ptr[i*(k-1)+j-1]+size;
            A_value.insert(A_value.end(),ne_X_value.begin()+begini, ne_X_value.begin()+endi);
            A_value.insert(A_value.end(),X_value.begin()+begini, X_value.begin()+endi);
            A_value.insert(A_value.end(),X_value.begin()+begini, X_value.begin()+endi);
            A_value.insert(A_value.end(),ne_X_value.begin()+begini, ne_X_value.begin()+endi);
            A_value.push_back(-1);
            for (L row = begini; row < endi; ++row)
                A_col_idx.push_back(X_col_idx[row]+2*(classyi-1)*d);
            for (L row = begini; row < endi; ++row)
                A_col_idx.push_back(X_col_idx[row]+(2*classyi-1)*d);
            for (L row = begini; row < endi; ++row)
                A_col_idx.push_back(X_col_idx[row]+2*j*d);
            for (L row = begini; row < endi; ++row)
                A_col_idx.push_back(X_col_idx[row]+(2*j+1)*d);
            A_col_idx.push_back(2*k*d+i);
        }
    }

    L mi = m;
    L me = 0;
    L nb = n;
    L nf = 0;

    std::string result_path = root_directory+ "/" + experiment + "/SVMtoLP";
    Libsvm_writeMat(result_path+"/svmtoLP",mi,n,A_row_ptr,A_col_idx,A_value,b);
    Makefile(result_path + "/Aeq");
    Makefile(result_path + "/beq");
    writeMeta(result_path +"/meta",mi,me,nb,nf);
    VectorWrite(result_path + "/c",c);
    writeMat(result_path +"/A",mi,n,A_row_ptr,A_col_idx,A_value);
    writeMat(result_path +"/Aeq",n);
    VectorWrite(result_path + "/b",b);



};



#endif //SDNA_LPUSE_H
