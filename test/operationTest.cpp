//
// Created by juzix on 2021/6/9.
//

#include "catch_amalgamated.hpp"
#include <bits/stdc++.h>
#include <typedef.hpp>
#include <fieldElement.hpp>
#include <groupElement.hpp>
#include <timer.hpp>
#include <wnaf.hpp>
#include <gmp.h>
#include <pippenger.hpp>

using std::vector;
using std::cerr;
using std::endl;
using namespace hyrax_p224;

class fieldElementTest : public fieldElement {
public:
    fieldElementTest(const mpz_t &other) {
        u32 dat[7];
        size_t countp = 7;
        mpz_export(dat, &countp, -1, sizeof (u32), 0, 0, other);

        u64 res = 0;
        u64 msk0 = (1 << 24) - 1;
        u64 msk1 = (1 << 16) - 1;
        u64 msk2 = (1 << 8) - 1;
        data[0] = (u64) dat[0] ^ ((u64) dat[1] & msk0) << 32;
        data[1] = ((u64) dat[1] >> 24) ^ ((u64) dat[2] << 8) ^ (((u64) dat[3] & msk1) << 40);
        data[2] = ((u64) dat[3] >> 16) ^ ((u64) dat[4] << 16) ^ (((u64) dat[5] & msk2) << 48);
        data[3] = ((u64) dat[5] >> 8) ^ ((u64) dat[6] << 24);
    }
};

typedef fieldElementTest F_2;

u32 mod_word_array[] = {0xffffffff, 0xffffffff, 0xffffffff, 0xffffffff, 0x00000000, 0x00000000, 0x00000001};
TEST_CASE("scalar_add_correctness", "[correctness]") {
    mpz_t md;
    mpz_init(md);
    mpz_import(md, 7, 1, sizeof (u32), 0, 0, mod_word_array);

    mpz_t a, b, c;
    mpz_init(a);
    mpz_init(b);
    mpz_init(c);

    int T = 100;
    gmp_randstate_t state;
    gmp_randinit_mt(state);

    while (T--) {
        mpz_urandomm(a, state, md);
        mpz_mod(a, a, md);
        mpz_urandomm(b, state, md);
        mpz_mod(b, b, md);
        mpz_add(c, a, b);
        mpz_mod(c, c, md);

//        gmp_printf ("a = %#56Zx\n", a);
//        gmp_printf ("b = %#56Zx\n", b);
//        gmp_printf ("c = %#56Zx\n", c);
//        gmp_printf ("md = %#56Zx\n", md);

        F_2 p224_a(a), p224_b(b), p224_c(c);
//        p224_a.print();
//        p224_b.print();
//        p224_c.print();
//        (p224_a + p224_b).print();
        REQUIRE((p224_a + p224_b) == p224_c);
    }

    mpz_clear(a);
    mpz_clear(b);
    mpz_clear(c);
    mpz_clear(md);
}

TEST_CASE("scalar_sub_correctness", "[correctness]") {
    mpz_t md;
    mpz_init(md);
    mpz_import(md, 7, 1, sizeof (u32), 0, 0, mod_word_array);

    mpz_t a, b, c;
    mpz_init(a);
    mpz_init(b);
    mpz_init(c);

    int T = 100;
    gmp_randstate_t state;
    gmp_randinit_mt(state);

    while (T--) {
        mpz_urandomm(a, state, md);
        mpz_mod(a, a, md);
        mpz_urandomm(b, state, md);
        mpz_mod(b, b, md);
        mpz_sub(c, a, b);
        mpz_mod(c, c, md);

//        gmp_printf ("a = %#56Zx\n", a);
//        gmp_printf ("b = %#56Zx\n", b);
//        gmp_printf ("c = %#56Zx\n", c);
//        gmp_printf ("md = %#56Zx\n", md);

        F_2 p224_a(a), p224_b(b), p224_c(c);
//        p224_a.print();
//        p224_b.print();
//        p224_c.print();
//        (p224_a - p224_b).print();
        REQUIRE(p224_a - p224_b == p224_c);
    }

    mpz_clear(a);
    mpz_clear(b);
    mpz_clear(c);
    mpz_clear(md);
}

TEST_CASE("scalar_mul_correctness", "[correctness]") {
    mpz_t md;
    mpz_init(md);
    mpz_import(md, 7, 1, sizeof (u32), 0, 0, mod_word_array);

    mpz_t a, b, c;
    mpz_init(a);
    mpz_init(b);
    mpz_init(c);

    int T = 10000;
    gmp_randstate_t state;
    gmp_randinit_mt(state);

    while (T--) {
        mpz_urandomm(a, state, md);
        mpz_mod(a, a, md);
        mpz_urandomm(b, state, md);
        mpz_mod(b, b, md);
        mpz_mul(c, a, b);
        mpz_mod(c, c, md);

        F_2 p224_a(a), p224_b(b), p224_c(c);
        REQUIRE(p224_a * p224_b == p224_c);
    }

    mpz_clear(a);
    mpz_clear(b);
    mpz_clear(c);
    mpz_clear(md);
}

TEST_CASE("exponential_test", "[correctness]") {
    groupElement bas = groupElement::random();
    fieldElement exp = fieldElement::random();

    REQUIRE(bas.isOnCurve());

    bool at_infinity = true;
    groupElement ans;
    for (u8 i = 223; i < 224; --i) {
        if (!at_infinity)
            ans.setDbl();

        if (exp.getBit(i)) {
            if (at_infinity) ans = bas;
            else ans += bas;
            at_infinity = false;
        }
    }

    if (at_infinity) ans.setInfinity();
    groupElement res = bas * exp;
    REQUIRE(res.isOnCurve());
    REQUIRE(ans == res);
}

TEST_CASE( "operation_test", "[correctness]" ) {
    const int N = 10;
    vector<fieldElement> sc_list(1 << N);
    for (auto &x: sc_list) x = fieldElement::random();

    vector<groupElement> pt_list(1 << N);
    for (auto &x: pt_list) x = groupElement::random();

    pippenger::mulExp mult(pt_list);
    auto res1 = mult.multiExponential(sc_list);
    auto res2 = pt_list[0] * sc_list[0];
    for (int j = 1; j < (1 << N); ++j)
        res2 += pt_list[j] * sc_list[j];

    REQUIRE(res1 == res2);
}

TEST_CASE("mul_time_test", "[time]") {
    int T = 10000;

    timer tm;
    while (T--) {
        fieldElement a = fieldElement::random(), b = fieldElement::random();
        tm.start();
        fieldElement c = a * b;
        tm.stop();
    }
    cerr << tm.elapse_sec() << endl;
}

TEST_CASE("dbl_time_test", "[time]") {
    int T = 10000;

    timer tm;
    while (T--) {
        groupElement g = groupElement::random();
        tm.start();
        g.setDbl();
        tm.stop();
    }
    cerr << tm.elapse_sec() << endl;
}

TEST_CASE( "operation_time_test", "[time]" ) {
    for (int i = 7; i <= 15; ++i) {
        vector<fieldElement> sc_list(1 << i);
        for (auto &x: sc_list) x = fieldElement::random();

        vector<groupElement> pt_list(1 << i);
        for (auto &x: pt_list) x = groupElement::random();

        pippenger::mulExp mult(pt_list);

        timer tm;
        tm.start();
        auto res1 = mult.multiExponential(sc_list);
        tm.stop();
//        auto res2 = pt_list[0] * sc_list[0];
//        for (int j = 0; j < (1 << i); ++j)
//            res2 += pt_list[j] * sc_list[j];
//        res1.print();
//        res2.print();
//        REQUIRE(res1 == res2);
        cerr << i << ", " << tm.elapse_sec() << endl;
    }
}

