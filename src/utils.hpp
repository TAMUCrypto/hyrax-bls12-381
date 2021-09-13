//
// Created by juzix on 2021/6/1.
//

#ifndef HYRAX_P224_UTILS_HPP
#define HYRAX_P224_UTILS_HPP

#include "typedef.hpp"
#include "fieldElement.hpp"

namespace hyrax_p224 {
    u8 myLog2(u64 x);

    bool checkPow2(u64 x);

    void split(vector<fieldElement> &L, vector<fieldElement> &R, const vector<fieldElement> &r);

    vector<fieldElement> expand(const vector<fieldElement> &v);

}



#endif //HYRAX_P224_UTILS_HPP
