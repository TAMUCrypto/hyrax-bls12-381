//
// Created by juzix on 2021/6/4.
//

#ifndef HYRAX_P224_POLYVERIFIER_HPP
#define HYRAX_P224_POLYVERIFIER_HPP


#include "polyProver.hpp"
#include "fieldElement.hpp"
#include "groupElement.hpp"
#include "utils.hpp"

namespace hyrax_p224 {
    using std::unique_ptr;
    class polyVerifier {
    public:
        polyVerifier(polyProver &_p, const vector<groupElement> &_gens);

        bool verify(const vector<fieldElement> &_x, const fieldElement &LZR);

        double getVT() { return vt.elapse_sec() - p.getPT(); }
    private:

        bool bulletVerify(vector<groupElement> g, vector<fieldElement> t, groupElement comm, fieldElement y);

        polyProver &p;

        vector<fieldElement> x, lx, rx;

        vector<groupElement> gens;
        groupElement blind;

        vector<groupElement> comm_Z;
        groupElement comm_ZR;

        timer vt;
        unique_ptr<pippenger::mulExp> commZ_multiplier;
    };
}


#endif //HYRAX_P224_POLYVERIFIER_HPP
