//
// Created by juzix on 2021/6/1.
//

#include "utils.hpp"
#include <vector>
#include <cmath>
using std::vector;

namespace hyrax_bls12_381 {
    u8 myLog2(u64 x) {
        return (u8) ceil(log(x) / log(2));
    }

    bool checkPow2(u64 x) {
        u8 tmp = (u8) round(log(x) / log(2));
        return x == 1ULL << tmp;
    }

    void split(vector<Fr> &L, vector<Fr> &R, const vector<Fr> &r) {
        u64 rsize = r.size() >> 1;
        u64 lsize = r.size() - rsize;
        L.clear();
        L.insert(L.end(), r.begin(), r.begin() + lsize);
        R.clear();
        R.insert(R.end(), r.begin() + lsize, r.end());
    }

    vector<Fr> expand(const vector<Fr> &v) {
        vector<Fr> V;
        vector<Fr> beta_f, beta_s;

        u64 first_half = v.size() >> 1;
        u64 second_half = v.size() - first_half;
        u64 mask = (1ULL << first_half) - 1;

        beta_f.resize(1ULL << first_half);
        beta_f[0] = Fr::one();
        for (u64 i = 0; i < first_half; ++i) {
            for (u64 j = 0; j < (1ULL << i); ++j) {
                auto tmp = beta_f.at(j) * v[i];
                beta_f[j | (1ULL << i)] = tmp;
                beta_f[j] = beta_f[j] - tmp;
            }
        }

        beta_s.resize(1ULL << second_half);
        beta_s[0] = Fr::one();
        for (u64 i = 0; i < second_half; ++i) {
            for (u64 j = 0; j < (1ULL << i); ++j) {
                auto tmp = beta_s[j] * v[(i + first_half)];
                beta_s[j | (1ULL << i)] = tmp;
                beta_s[j] = beta_s[j] - tmp;
            }
        }

        u64 size = 1ULL << v.size();
        V.resize(size);
        for (u64 i = 0; i < size; ++i)
            V[i] = beta_f[i & mask] * beta_s[i >> first_half];
        return V;
    }
}
