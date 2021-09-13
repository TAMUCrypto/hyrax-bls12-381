//
// Created by juzix on 2021/6/11.
//

#ifndef HYRAX_P224_WNAF_HPP
#define HYRAX_P224_WNAF_HPP

#include <bits/stdc++.h>
#include "typedef.hpp"
#include "groupElement.hpp"

using std::vector;

namespace hyrax_p224 {
namespace wnaf {
    void ec_compute_wNAF(i8 *out, const fieldElement &scalar, size_t bits, int w);

    void compute_precomp(groupElement *out, const groupElement p, size_t len);

    void lookup_precomp(groupElement &out, const groupElement *precomp, int digit);

    void mul_batch(groupElement &r, const fieldElement *g_scalar, const groupElement *points, const fieldElement *scalars, size_t num);

    groupElement vartimeMultiscalarMul(const vector<fieldElement> &exp, const vector<groupElement> &bas);

    groupElement
    vartimeMultiscalarMul(const vector<fieldElement>::const_iterator &exp, const vector<groupElement>::const_iterator &bas, u64 n);
}

}
#endif //HYRAX_P224_WNAF_HPP
