//
// Created by juzix on 2021/6/4.
//

#include "polyVerifier.hpp"

namespace hyrax_p224 {
    using std::unique_ptr;
    using std::make_unique;
    polyVerifier::polyVerifier(polyProver &_p, const vector<groupElement> &_gens) : p(_p), gens(_gens) {
        timer tmp_timer;
        tmp_timer.start();
        vt.start();
        comm_Z = p.commit();
        commZ_multiplier = make_unique<pippenger::mulExp>(comm_Z);
        vt.stop();
        tmp_timer.stop();
        fprintf(stderr, "commit time: %.4f\n", tmp_timer.elapse_sec());
    }

    bool polyVerifier::verify(const vector<fieldElement> &_x, const fieldElement &LZR) {
        fprintf(stderr, "Poly commit for 2^%d input.\n", (int) _x.size());
        vt.start();
        x = _x;
        split(lx, rx, x);
        p.initBulletProve(lx, rx);

        comm_ZR = commZ_multiplier->multiExponential(expand(rx));

        bool res = bulletVerify(gens, lx, comm_ZR, LZR);
        vt.stop();

        return res;
    }

    bool polyVerifier::bulletVerify(vector<groupElement> g, vector<fieldElement> t, groupElement comm,
                                    fieldElement y) {
        timer tmp_timer;
        tmp_timer.start();
        groupElement lcomm, rcomm;
        fieldElement ly, ry;

        assert(checkPow2(g.size()));
        auto logn = t.size();
        while (true) {
            p.bulletProve(lcomm, rcomm, ly, ry);

            fieldElement randomness = fieldElement::random();
            fieldElement irandomness = randomness.inv();

            p.bulletUpdate(randomness);

            u64 hsize = g.size() >> 1;
            for (u64 i = 0; i < hsize; ++i)
                g[i] = g[i] * irandomness + g[i + hsize];
            g.resize(hsize);

            comm = lcomm * randomness + comm + rcomm * irandomness;
            if (y != ly * (fieldElement::one() - t.back()) + ry * t.back()) {
                fprintf(stderr, "y incorrect at %d.\n", (int) (logn - t.size()));
                return false;
            }
            y = ly * randomness + ry;

            if (t.size() == 1) {
                bool res = p.bulletOpen() == y;

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