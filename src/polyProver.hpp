//
// Created by juzix on 2021/6/4.
//

#ifndef HYRAX_P224_POLYPROVER_HPP
#define HYRAX_P224_POLYPROVER_HPP


#include "timer.hpp"
#include "utils.hpp"
#include "typedef.hpp"
#include <vector>
#include <mcl/bls12_381.hpp>
using std::vector;
using namespace mcl::bn;

namespace hyrax_bls12_381 {
    class polyProver {
    public:
        polyProver(const vector<Fr> &_Z, const vector<G1> &_gens);

        vector<G1> commit();

        Fr evaluate(const vector<Fr> &x);

        double getPT() const;

        double getPS() const;

        void initBulletProve(const vector<Fr> &_lx, const vector<Fr> &_rx);

        void bulletProve(G1 &lcomm, G1 &rcomm, Fr &ly, Fr &ry);

        void bulletUpdate(const Fr &randomness);

        Fr bulletOpen();

    private:
        vector<G1> gens;
    public:
        const vector<G1> &getGens() const;

    private:
        vector<G1> comm_Z;

        vector<Fr> Z, L, R, t;
        vector<Fr> RZ;
        Fr scale;
        u8 bit_length;
        timer pt;   // s
        u64 ps;     // KB

        vector<G1> bullet_g;
        vector<Fr> bullet_a;
    };
}


#endif //HYRAX_P224_POLYPROVER_HPP
