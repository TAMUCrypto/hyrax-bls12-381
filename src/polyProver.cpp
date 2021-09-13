//
// Created by juzix on 2021/6/4.
//

#include "polyProver.hpp"
#include "wnaf.hpp"

using std::vector;
namespace hyrax_p224 {
    polyProver::polyProver(const vector<fieldElement> &_Z, const vector<groupElement> &_gens) :
            Z(_Z), gens(_gens), multiplier(_gens) {
        bit_length = myLog2(_Z.size());
        ps = 0;
    }

    vector<groupElement> polyProver::commit() {
        pt.start();
        u8 r_bit_length = bit_length >> 1;
        u8 l_bit_length = bit_length - r_bit_length;

        u64 rsize = 1ULL << r_bit_length, lsize = 1ULL << l_bit_length;
        assert(lsize == gens.size());

        comm_Z.resize(rsize);
        for (u64 i = 0; i < rsize; ++i)
            comm_Z[i] = multiplier.multiExponential(Z.begin() + i * lsize, Z.begin() + (i + 1) * lsize);

        pt.stop();
        ps += (3 * sizeof(fieldElement)) * comm_Z.size();
        return comm_Z;
    }

    fieldElement polyProver::evaluate(const vector<fieldElement> &x) {
        auto X = expand(x);
        auto res = fieldElement::zero();
        assert(X.size() == Z.size());
        for (u64 i = 0; i < X.size(); ++i) res += Z[i] * X[i];
        return res;
    }

    double polyProver::getPT() const {
        return pt.elapse_sec();
    }

    double polyProver::getPS() const {      // KB
        return ps / 1024.0;
    }

    void polyProver::initBulletProve(const vector<fieldElement> &_lx, const vector<fieldElement> &_rx) {
        pt.start();

        t = _lx;
        L = expand(_lx);
        R = expand(_rx);

        u64 lsize_ex = L.size(), rsize_ex = R.size();
        assert(lsize_ex * rsize_ex == Z.size());
        assert(lsize_ex == gens.size());

        ZR.resize(lsize_ex, fieldElement::zero());
        for (u64 j = 0; j < rsize_ex; ++j)
            for (u64 i = 0; i < lsize_ex; ++i)
                ZR[i] += Z[j * lsize_ex + i] * R[j];

//        pippenger::mulExp m2(comm_Z);
//        assert(multiplier.multiExponential(ZR) == m2.multiExponential(R));
        groupElement res;
        for (int i = 0; i < R.size(); ++i) {
            for (int j = 0; j < L.size(); ++j) {
                if (!i && !j) res = gens[j] * (R[i] * Z[i * lsize_ex + j]);
                else res += gens[j] * (R[i] * Z[i * lsize_ex + j]);
            }
        }
        auto res2 = comm_Z[0] * R[0];
        for (int i = 1; i < R.size(); ++i) res2 += comm_Z[i] * R[i];
        res.print();
        res2.print();
        assert(res2 == res);

        bullet_g = gens;
        bullet_a = ZR;
        scale = fieldElement::one();
        pt.stop();
    }

    void polyProver::bulletProve(groupElement &lcomm, groupElement &rcomm, fieldElement &ly, fieldElement &ry) {
        pt.start();
        assert(!(bullet_a.size() & 1));
        u64 hsize = bullet_a.size() >> 1;

        pippenger::mulExp bullet_multiplier1(vector<groupElement>(bullet_g.begin(), bullet_g.begin() + hsize));
        pippenger::mulExp bullet_multiplier2(vector<groupElement>(bullet_g.begin() + hsize, bullet_g.end()));
        lcomm = bullet_multiplier2.multiExponential(bullet_a.begin(), bullet_a.begin() + hsize);
        rcomm = bullet_multiplier1.multiExponential(bullet_a.begin() + hsize, bullet_a.end());

        scale *= (fieldElement::one() - t.back()).inv();
        ly = fieldElement::innerProd(bullet_a.begin(), L.begin(), hsize) * scale;
        ry = fieldElement::innerProd(bullet_a.begin() + hsize, L.begin(), hsize) * scale;
        pt.stop();
        ps += ((3 * sizeof(fieldElement)) + sizeof(fieldElement)) * 2;
    }

    void polyProver::bulletUpdate(const fieldElement &randomness) {
        pt.start();
        auto irandomness = randomness.inv();
        u64 hsize = bullet_a.size() >> 1;
        for (u64 i = 0; i < hsize; ++i) bullet_a[i] = bullet_a[i] * randomness + bullet_a[i + hsize];
        for (u64 i = 0; i < hsize; ++i) bullet_g[i] = bullet_g[i] * irandomness + bullet_g[i + hsize];
        bullet_a.resize(hsize);
        bullet_g.resize(hsize);
        t.pop_back();
        pt.stop();
    }

    fieldElement polyProver::bulletOpen() {
        assert(bullet_a.size() == 1);

        ps += sizeof(fieldElement);
        return bullet_a.back();
    }

    const vector<groupElement> &polyProver::getGens() const {
        return gens;
    }
}