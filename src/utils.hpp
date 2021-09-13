//
// Created by juzix on 2021/6/1.
//

#ifndef HYRAX_P224_UTILS_HPP
#define HYRAX_P224_UTILS_HPP

#include "typedef.hpp"
#include <vector>
#include <mcl/bls12_381.hpp>

using std::vector;
using namespace mcl::bn;
namespace hyrax_bls12_381 {
    u8 myLog2(u64 x);

    bool checkPow2(u64 x);

    void split(vector<Fr> &L, vector<Fr> &R, const vector<Fr> &r);

    vector<Fr> expand(const vector<Fr> &v);

}



#endif //HYRAX_P224_UTILS_HPP
