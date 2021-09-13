//
// Created by juzix on 2021/6/17.
//

#include "typedef.hpp"
#include "fieldElement.hpp"
#include "polyProver.hpp"
#include "polyVerifier.hpp"
#include <bits/stdc++.h>

using std::vector;
using std::string;
using namespace hyrax_p224;

#define LOGN_ID 0
#define PT_ID 1
#define VT_ID 2
#define PS_ID 3

vector<string> ans(4);

template <typename T>
string to_string_wp(const T a_value, const int n = 4) {
    std::ostringstream out;
    out.precision(n);
    out << std::fixed << a_value;
    return out.str();
}

int main(int argc, char *argv[]) {
    fieldElement::init();

    // Evaluate the multiplication
    timer mul_timer;
    int mul_repno = 10000;
    for (int i = 0; i < mul_repno; ++i) {
        fieldElement a, b, c;
        a = fieldElement::random();
        b = fieldElement::random();
        mul_timer.start();
        c = a * b;
        mul_timer.stop();
        if (c != (a * b)) puts("Wrong!");
    }
    printf("mul time for %d opts: %.5f(s)\n", mul_repno, mul_timer.elapse_sec());

    // Evaluate the exponential
    printf("Evaluate the exponential\n");
    printf("logn,\ttime\n");
    int exp_repno = 100;
    for (int i = 7; i < 13; ++i) {
        timer exp_timer;
        vector<fieldElement> a(1 << i);
        vector<groupElement> b(1 << i);
        groupElement c;
        for (int j = 0; j < exp_repno; ++j) {
            for (auto &x: a) x = fieldElement::random();
            for (auto &x: b) x = groupElement::random();
            exp_timer.start();
            pippenger::mulExp multiplier(b);
            c = multiplier.multiExponential(a.begin(), a.end());
            exp_timer.stop();
        }
        printf("%d,\t%.5f\n", i, exp_timer.elapse_sec() / exp_repno);
    }

    // Evaluate the polynomial
    printf("Evaluate the polynomial\n");
    printf("logn,\tpt,\tvt,\tps\n");
    for (int logn = 15; logn < 22; ++logn) {
        u64 n = 1ULL << logn;
        u64 n_sqrt = 1ULL << (logn - (logn >> 1));
        vector<fieldElement> poly_coeff(n);
        vector<groupElement> gens(n_sqrt);
        vector<fieldElement> r(logn);

        for (u64 i = 0; i < n; ++i) poly_coeff[i] = fieldElement::random();
        for (u64 i = 0; i < n_sqrt; ++i) gens[i] = groupElement::random();
        for (u8 i = 0; i < logn; ++i) r[i] = fieldElement();
        hyrax_p224::polyProver p(poly_coeff, gens);
        hyrax_p224::polyVerifier v(p, gens);
        if (!v.verify(r, p.evaluate(r))) puts("Wrong!");

        ans[LOGN_ID] = std::to_string(logn);
        ans[PT_ID] = to_string_wp(p.getPT());
        ans[VT_ID] = to_string_wp(v.getVT());
        ans[PS_ID] = to_string_wp(p.getPS());

        for (auto & an : ans)
            printf("%s,\t", an.c_str());
        puts("");
    }
    return 0;
}
