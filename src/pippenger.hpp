//
// Created by juzix on 2021/6/16.
//

#ifndef HYRAX_P224_PIPPENGER_HPP
#define HYRAX_P224_PIPPENGER_HPP


#include "groupElement.hpp"

namespace hyrax_p224 {
namespace pippenger {
    struct mulExp {
    public:
        mulExp(const vector<groupElement> &_gens);
        mulExp(const vector<groupElement> &_gens, u64 _tabl_size, u64 _win_size);

        groupElement multiExponential(const vector<fieldElement> &exp);
        groupElement multiExponential(const vector<fieldElement>::const_iterator &first,
                                      const vector<fieldElement>::const_iterator &last);

        groupElement multiExponential(const vector<fieldElement>::const_iterator &first,
                                      const vector<fieldElement>::const_iterator &last,
                                      u64 gen_first, u64 gen_last);

    private:
        vector<groupElement> gens;
        vector<groupElement> gens_ext; // g_i^{2^j}
        vector<vector<groupElement>> tabl;

        const u64 N;
        static const u64 lambda = 224;
        u64 TABL_SIZE;
        u64 TABL_CNT;
        u64 WIN_SIZE;

    };
}
}

#endif //HYRAX_P224_PIPPENGER_HPP
