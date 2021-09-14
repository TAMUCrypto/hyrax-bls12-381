//
// Created by juzix on 2021/6/4.
//

#include "polyProver.hpp"
#include <mcl/bls12_381.hpp>

#define G1_SIZE (Fp::getByteSize())
#define Fr_SIZE (Fr::getByteSize())

using std::vector;
namespace hyrax_bls12_381 {
    polyProver::polyProver(const vector<Fr> &_Z, const vector<G1> &_gens) :
            Z(_Z), gens(_gens) {
        bit_length = myLog2(_Z.size());
        ps = 0;
    }

    vector<G1> polyProver::commit() {
        pt.start();
        u8 r_bit_length = bit_length >> 1;
        u8 l_bit_length = bit_length - r_bit_length;

        u64 rsize = 1ULL << r_bit_length, lsize = 1ULL << l_bit_length;
        assert(lsize == gens.size());

        comm_Z.resize(rsize);
        for (u64 i = 0; i < rsize; ++i)
            G1::mulVec(comm_Z[i], gens.data(), Z.data() + i * lsize, lsize);

        pt.stop();
        ps += G1_SIZE * comm_Z.size();
        return comm_Z;
    }

    Fr polyProver::evaluate(const vector<Fr> &x) {
        auto X = expand(x);
        Fr res;
        assert(X.size() == Z.size());
        for (u64 i = 0; i < X.size(); ++i) res = !i ? Z[i] * X[i] : res + Z[i] * X[i];
        return res;
    }

    double polyProver::getPT() const {
        return pt.elapse_sec();
    }

    double polyProver::getPS() const {      // KB
        return ps / 1024.0;
    }

    void polyProver::initBulletProve(const vector<Fr> &_lx, const vector<Fr> &_rx) {
        Fr zero;
        zero.clear();
        pt.start();

        t = _lx;
        L = expand(_lx);
        R = expand(_rx);

        u64 lsize_ex = L.size(), rsize_ex = R.size();
        assert(lsize_ex * rsize_ex == Z.size());
        assert(lsize_ex == gens.size());

        RZ.resize(lsize_ex, zero);
        for (u64 j = 0; j < rsize_ex; ++j)
            for (u64 i = 0; i < lsize_ex; ++i)
                RZ[i] = !j ? R[j] * Z[j * lsize_ex + i] : RZ[i] + R[j] * Z[j * lsize_ex + i];

        bullet_g = gens;
        bullet_a = RZ;
        scale = Fr::one();
        pt.stop();
    }

    void polyProver::bulletProve(G1 &lcomm, G1 &rcomm, Fr &ly, Fr &ry) {
        pt.start();
        assert(!(bullet_a.size() & 1));
        u64 hsize = bullet_a.size() >> 1;

        G1::mulVec(lcomm, bullet_g.data(), bullet_a.data(), hsize);
        G1::mulVec(rcomm, bullet_g.data() + hsize, bullet_a.data() + hsize, hsize);

        Fr tmp;
        Fr::inv(tmp, Fr::one() - t.back());
        scale *= tmp;

        for (int i = 0; i < hsize; ++i) {
            ly = !i ? bullet_a[i] * L[i] : ly + bullet_a[i] * L[i];
            ry = !i ? bullet_a[i + hsize] * L[i] : ry + bullet_a[i + hsize] * L[i];
        }
        ly *= scale;
        ry *= scale;
        pt.stop();
        ps += (G1_SIZE + Fr_SIZE) * 2;
    }

    void polyProver::bulletUpdate(const Fr &randomness) {
        pt.start();
        Fr irandomness;
        Fr::inv(irandomness, randomness);
        u64 hsize = bullet_a.size() >> 1;
        for (u64 i = 0; i < hsize; ++i) bullet_a[i] = bullet_a[i] * randomness + bullet_a[i + hsize];
        for (u64 i = 0; i < hsize; ++i) bullet_g[i] = bullet_g[i] * irandomness + bullet_g[i + hsize];
        bullet_a.resize(hsize);
        bullet_g.resize(hsize);
        t.pop_back();
        pt.stop();
    }

    Fr polyProver::bulletOpen() {
        assert(bullet_a.size() == 1);

        ps += Fr_SIZE;
        return bullet_a.back();
    }

    const vector<G1> &polyProver::getGens() const {
        return gens;
    }
}