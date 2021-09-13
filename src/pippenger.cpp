//
// Created by juzix on 2021/6/16.
//

#include "pippenger.hpp"

namespace hyrax_p224 {
namespace pippenger {

    mulExp::mulExp(const vector<groupElement> &_gens): gens(_gens), N(_gens.size()) {
        TABL_SIZE = 1;
        TABL_CNT = TABL_SIZE * _gens.size();
        WIN_SIZE = std::max(6, int( (log(N) - log(log(N) / log(2.))) / log(2.) ));
        for (u64 i = 0; i < N; ++i) {
            gens_ext.push_back(gens[i]);
            for (int j = 1; j < TABL_SIZE; ++j) gens_ext.push_back(gens_ext.back().dbl());
        }

        tabl.resize((gens_ext.size() + WIN_SIZE - 1) / WIN_SIZE);
        for (auto &x: tabl) x.resize(1ULL << WIN_SIZE);
        for (u64 i = 0; i < gens_ext.size(); ++i) {
            u64 x = i / WIN_SIZE;
            u64 y = i % WIN_SIZE;

            u64 pos = (1ULL << y);
            tabl[x][pos] = gens_ext[i];
            for (u64 j = 1; j < pos; ++j)
                tabl[x][j ^ pos] = tabl[x][j] + gens_ext[i];
        }
    }

    groupElement mulExp::multiExponential(const vector<fieldElement> &exp) {
        return multiExponential(exp.begin(), exp.end(), 0, gens.size());
    }

    groupElement mulExp::multiExponential(const vector<fieldElement>::const_iterator &first,
                                          const vector<fieldElement>::const_iterator &last) {
        return multiExponential(first, last, 0, gens.size());
    }

    groupElement mulExp::multiExponential(const vector<fieldElement>::const_iterator &first,
                                          const vector<fieldElement>::const_iterator &last,
                                          u64 gen_first, u64 gen_last) {
        assert(last - first == N);

        vector<groupElement> G;
        u8 max_bit = 0;
        for (u64 i = gen_first; i < gen_last; ++i) {
            max_bit = std::max(max_bit, (/*bit_widths[i] = */first[i].getBitWidth()));
        }

        if (!max_bit) return groupElement() * fieldElement::zero();
        G.resize((max_bit + TABL_SIZE - 1) / TABL_SIZE);
        for (u8 i = 0; i < max_bit; i += TABL_SIZE) {
            bool at_infinity = true;
            u64 cnt = (gen_first * TABL_SIZE) % WIN_SIZE;
            u64 x = (gen_first * TABL_SIZE) / WIN_SIZE, y = 0;
            for (u64 j = gen_first; j < gen_last; ++j)
                    for (u64 k = 0; k < TABL_SIZE; ++k) {
                        y = y ^ (first[j - gen_first].getBit(i + k) << cnt);
                        ++cnt;
                        if (cnt == WIN_SIZE || j == gen_last - 1 && k == TABL_SIZE - 1) {
                            if (y) {
                                if (at_infinity) G[i / TABL_SIZE] = tabl[x][y];
                                else G[i / TABL_SIZE] += tabl[x][y];
                                at_infinity = false;
                            }
                            cnt = 0;
                            y = 0;
                            ++x;
                        }
                    }
        }

        groupElement res = G.back();
        for (u64 i = G.size() - 2; i < G.size() - 1; --i) {
            for (u8 j = 0; j < TABL_SIZE; ++j) res.setDbl();
            res += G[i];
        }

        return res;
    }

}
}