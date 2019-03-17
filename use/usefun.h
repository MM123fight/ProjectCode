//
// Created by Lu Meng on 2018/10/23.
//

#ifndef USE_FUN_H
#define USE_FUN_H

#include "stdafx.h"
//Transfer number to the type of char or string
template <typename T >
std::string to_string(const T t){
    std::stringstream s;
    s << t;
    return s.str();
}
#endif //USE_FUN_H
