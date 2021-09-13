//
// Created by juzix on 2021/6/4.
//

#ifndef HYRAX_P224_POLYVERIFIER_HPP
#define HYRAX_P224_POLYVERIFIER_HPP


#include "polyProver.hpp"
#include "utils.hpp"
#include <vector>
using std::vector;

namespace hyrax_bls12_381 {
    class polyVerifier {
    public:
        polyVerifier(polyProver &_p, const vector<G1> &_gens);

        bool verify(const vector<Fr> &_x, const Fr &RZL);

        double getVT() { return vt.elapse_sec() - p.getPT(); }
    private:

        bool bulletVerify(vector<G1> g, vector<Fr> t, G1 comm, Fr y);

        polyProver &p;

        vector<Fr> x, lx, rx;

        vector<G1> gens;
        G1 blind;

        vector<G1> comm_Z;
        G1 comm_RZ;

        timer vt;
    };
}


#endif //HYRAX_P224_POLYVERIFIER_HPP
