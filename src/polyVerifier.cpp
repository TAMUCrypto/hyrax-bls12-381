//
// Created by juzix on 2021/6/4.
//

#include "polyVerifier.hpp"

namespace hyrax_bls12_381 {
    polyVerifier::polyVerifier(polyProver &_p, const vector<G1> &_gens) : p(_p), gens(_gens) {
        timer tmp_timer;
        tmp_timer.start();
        vt.start();
        comm_Z = p.commit();
        vt.stop();
        tmp_timer.stop();
        fprintf(stderr, "commit time: %.4f\n", tmp_timer.elapse_sec());
    }

    bool polyVerifier::verify(const vector<Fr> &_x, const Fr &RZL) {
        fprintf(stderr, "Poly commit for 2^%d input.\n", (int) _x.size());
        vt.start();
        x = _x;
        split(lx, rx, x);
        p.initBulletProve(lx, rx);

        auto R = expand(rx);
        assert(comm_Z.size() == R.size());
        G1::mulVec(comm_RZ, comm_Z.data(), R.data(), comm_Z.size());
        bool res = bulletVerify(gens, lx, comm_RZ, RZL);
        vt.stop();

        return res;
    }

    bool polyVerifier::bulletVerify(vector<G1> g, vector<Fr> t, G1 comm,
                                    Fr y) {
        timer tmp_timer;
        tmp_timer.start();
        G1 lcomm, rcomm;
        Fr ly, ry;

        assert(checkPow2(g.size()));
        auto logn = t.size();
        while (true) {
            p.bulletProve(lcomm, rcomm, ly, ry);

            Fr randomness;
            randomness.setByCSPRNG();
            Fr irandomness;
            Fr::inv(irandomness, randomness);

            p.bulletUpdate(randomness);

            u64 hsize = g.size() >> 1;
            for (u64 i = 0; i < hsize; ++i)
                g[i] = g[i] * irandomness + g[i + hsize];
            g.resize(hsize);

            comm = lcomm * randomness + comm + rcomm * irandomness;
            if (y != ly * (Fr::one() - t.back()) + ry * t.back()) {
                fprintf(stderr, "y incorrect at %d.\n", (int) (logn - t.size()));
                return false;
            }
            y = ly * randomness + ry;

            if (t.size() == 1) {
                bool res = p.bulletOpen() == y && comm == g.back() * y;

                tmp_timer.stop();
                fprintf(stderr, "bulletProve time: %.4f\n", tmp_timer.elapse_sec());

                if (!res) {
                    fprintf(stderr, "last step incorrect.\n");
                    return false;
                }
                return true;
            }
            t.pop_back();
        }
    }
}