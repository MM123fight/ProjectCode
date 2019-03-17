//
// Created by Lu Meng on 2018/11/27.
//
#include "../Problem/ProblemHeader.h"
#include "../LossFun/LossFunHeader.h"
#include "AbsAlgorithm.h"
#include "test.h"
#include "Abstest.h"

#include "ParamAdapPPA.h"

#include <algorithm>
#include <iostream>
#include <vector>
#include <cmath>


static bool abs_compare(int a, int b)
{
    return (std::abs(a) < std::abs(b));
}

int main()
{
    double varepsilon;
    double varepsilon_initial = 1.;
    double inv_idx = 1.1;
    unsigned int iter =1;
    Print("check");
    Print("varepsilon", varepsilon);
    inv_update(varepsilon,varepsilon_initial,inv_idx,iter);
    Print("varepsilon", varepsilon);
    std::vector<int> x;
    x = std::vector<int>(8);
    Print(x);
    int a = ceil(11);
    std::cout << a << std::endl;
    std::vector<int> v{ 3, 1, -14, 1, 5, 9 };
    std::vector<int>::iterator result;

    A* inst = new B();
    inst->print();
    inst->print();

    result = std::max_element(v.begin(), v.end());
    std::cout << "max element at: " << std::distance(v.begin(), result) << '\n';

    result = std::max_element(v.begin(), v.end(), abs_compare);
    std::cout << "max element (absolute) at: " << std::distance(v.begin(), result);

    delete inst;

    return 0;
}
/*
int main() {
    int p1 =10;
    int *p2 = new int(8);
    int *p3(p2);
    p3 = p2;

    Print(*p2);
    Print(*p3);

    delete p3;
    Print(*p2);
    Print(*p3);
    //delete p2;







    return 0;
}
 */