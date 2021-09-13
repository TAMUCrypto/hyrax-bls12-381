//
// Created by juzix on 2021/6/8.
//

#include "catch_amalgamated.hpp"
#include <bits/stdc++.h>
#include <typedef.hpp>
#include <polyProver.hpp>
#include <polyVerifier.hpp>
#include <mcl/bls12_381.hpp>

#define CATCH_CONFIG_MAIN

using std::vector;
using namespace hyrax_bls12_381;
using namespace mcl::bn;

TEST_CASE( "poly_correctness", "[correctness]" ) {
    initPairing(mcl::BLS12_381);
    u8 logn = 15;
    u64 n = 1ULL << logn;
    u64 n_sqrt = 1ULL << (logn - (logn >> 1));
    vector<Fr> poly_coeff(n);
    vector<G1> gens(n_sqrt);
    vector<Fr> r(logn);

    for (u64 i = 0; i < n; ++i) poly_coeff[i].setByCSPRNG();
    for (u64 i = 0; i < n_sqrt; ++i) {
        Fr tmp;
        tmp.setByCSPRNG();
        gens[i] = mcl::bn::getG1basePoint() * tmp;
    }
    for (u8 i = 0; i < logn; ++i) r[i].setByCSPRNG();
    polyProver p(poly_coeff, gens);
    polyVerifier v(p, gens);
    REQUIRE(v.verify(r, p.evaluate(r)));
}

//TEST_CASE( "poly_time", "[time]" ) {
//    initPairing(mcl::BLS12_381);
//    for (u8 logn = 15; logn < 30; ++logn) {
//        u64 n = 1ULL << logn;
//        u64 n_sqrt = 1ULL << (logn - (logn >> 1));
//        vector<Fr> poly_coeff(n);
//        vector<G1> gens(n_sqrt);
//        vector<Fr> r(logn);
//
//        for (u64 i = 0; i < n; ++i) poly_coeff[i].setByCSPRNG();
//        for (u64 i = 0; i < n_sqrt; ++i) {
//            Fr tmp;
//            tmp.setByCSPRNG();
//            gens[i] = mcl::bn::getG1basePoint() * tmp;
//        }
//        for (u8 i = 0; i < logn; ++i) r[i].setByCSPRNG();
//        polyProver p(poly_coeff, gens);
//        polyVerifier v(p, gens);
//        REQUIRE(v.verify(r, p.evaluate(r)));
//        printf("logn = %d, pt = %.3f, vt = %.3f, ps = %.3f\n", (int) logn, p.getPT(), v.getVT(), p.getPS());
//    }
//}