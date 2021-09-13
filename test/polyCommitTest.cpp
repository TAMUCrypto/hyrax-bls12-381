//
// Created by juzix on 2021/6/8.
//

#include "catch_amalgamated.hpp"
#include <bits/stdc++.h>
#include <typedef.hpp>
#include <fieldElement.hpp>
#include <polyProver.hpp>
#include <polyVerifier.hpp>

#define CATCH_CONFIG_MAIN

using std::vector;
using namespace hyrax_p224;

TEST_CASE( "poly_correctness", "[correctness]" ) {
    u8 logn = 15;
    u64 n = 1ULL << logn;
    u64 n_sqrt = 1ULL << (logn - (logn >> 1));
    vector<fieldElement> poly_coeff(n);
    vector<groupElement> gens(n_sqrt);
    vector<fieldElement> r(logn);

    for (u64 i = 0; i < n; ++i) poly_coeff[i] = fieldElement::random();
    for (u64 i = 0; i < n_sqrt; ++i) gens[i] = groupElement::random();
    for (u8 i = 0; i < logn; ++i) r[i] = fieldElement::random();
    polyProver p(poly_coeff, gens);
    polyVerifier v(p, gens);
    REQUIRE(v.verify(r, p.evaluate(r)));
}

//TEST_CASE( "poly_time", "[time]" ) {
//    for (u8 logn = 15; logn < 30; ++logn) {
//        u64 n = 1ULL << logn;
//        u64 n_sqrt = 1ULL << (logn - (logn >> 1));
//        vector<fieldElement> poly_coeff(n);
//        vector<groupElement> gens(n_sqrt);
//        vector<fieldElement> r(logn);
//
//        for (u64 i = 0; i < n; ++i) poly_coeff[i] = fieldElement::random();
//        for (u64 i = 0; i < n_sqrt; ++i) gens[i] = groupElement::random();
//        for (u8 i = 0; i < logn; ++i) r[i] = fieldElement::random();
//        polyProver p(poly_coeff, gens);
//        polyVerifier v(p, gens);
//        REQUIRE(v.verify(r, p.evaluate(r)));
//        printf("logn = %d, pt = %.3f, vt = %.3f, ps = %.3f\n", (int) logn, p.getPT(), v.getVT(), p.getPS());
//    }
//}