//
// Created by juzix on 2021/6/4.
//

#ifndef HYRAX_P224_POLYPROVER_HPP
#define HYRAX_P224_POLYPROVER_HPP


#include "fieldElement.hpp"
#include "timer.hpp"
#include "groupElement.hpp"
#include "utils.hpp"
#include "pippenger.hpp"

namespace hyrax_p224 {
    class polyProver {
    public:
        polyProver(const vector<fieldElement> &_Z, const vector<groupElement> &_gens);

        vector<groupElement> commit();

        fieldElement evaluate(const vector<fieldElement> &x);

        double getPT() const;

        double getPS() const;

        void initBulletProve(const vector<fieldElement> &_lx, const vector<fieldElement> &_rx);

        void bulletProve(groupElement &lcomm, groupElement &rcomm, fieldElement &ly, fieldElement &ry);

        void bulletUpdate(const fieldElement &randomness);

        fieldElement bulletOpen();

    private:
        vector<groupElement> gens;
    public:
        const vector<groupElement> &getGens() const;

    private:
        vector<groupElement> comm_Z;

        vector<fieldElement> Z, L, R, t;
        vector<fieldElement> ZR;
        fieldElement scale;
        u8 bit_length;
        timer pt;   // s
        u64 ps;     // KB

        vector<groupElement> bullet_g;
        vector<fieldElement> bullet_a;

        pippenger::mulExp multiplier;
    };
}


#endif //HYRAX_P224_POLYPROVER_HPP
