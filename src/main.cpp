//
// Created by juzix on 2021/6/17.
//

#include "typedef.hpp"
#include "polyProver.hpp"
#include "polyVerifier.hpp"
#include <bits/stdc++.h>
#include <mcl/bls12_381.hpp>

using std::vector;
using std::string;
using namespace hyrax_bls12_381;
using namespace mcl::bn;

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
    u8 logn = argc == 1 ? 10 : atoi(argv[1]);

    initPairing(mcl::BLS12_381);
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
    hyrax_bls12_381::polyProver p(poly_coeff, gens);
    hyrax_bls12_381::polyVerifier v(p, gens);
    assert(v.verify(r, p.evaluate(r)));

    ans[LOGN_ID] = std::to_string(logn);
    ans[PT_ID] = to_string_wp(p.getPT());
    ans[VT_ID] = to_string_wp(v.getVT());
    ans[PS_ID] = to_string_wp(p.getPS());

    for (int i = 0; i < ans.size(); ++i)
        printf("%s, ", ans[i].c_str());
    puts("");
    return 0;
}
