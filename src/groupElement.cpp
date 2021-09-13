//
// Created by juzix on 2021/6/1.
//

#include "groupElement.hpp"
#include "fieldElement.hpp"
#include "p224.hpp"

using std::move;

namespace hyrax_p224 {
groupElement groupElement::random() {
    return groupElement(__gens) * fieldElement::random();
}

groupElement groupElement::zero() {
    groupElement res;

    return res;
}

groupElement &groupElement::operator=(const groupElement &other) {
    p224_felem_assign(x, other.x);
    p224_felem_assign(y, other.y);
    p224_felem_assign(z, other.z);
    return *this;
}

groupElement groupElement::operator+(const groupElement &other) const {
    groupElement res;
    p224_point_add(res.x, res.y, res.z, x, y, z, 0, other.x, other.y, other.z);
    p224_felem_contract(res.x, res.x);
    p224_felem_contract(res.y, res.y);
    p224_felem_contract(res.z, res.z);
    return res;
}

groupElement &groupElement::operator+=(const groupElement &other) {
    *this = *this + other;
    return *this;
}

void ec_GFp_nistp_recode_scalar_bits(u32 *sign, u32 *digit,
                                     u32 in) {
    u32 s, d;

    s = ~((in >> 5) - 1); /* sets all bits to MSB(in), 'in' seen as
                          * 6-bit value */
    d = (1 << 6) - in - 1;
    d = (d & s) | (in & ~s);
    d = (d >> 1) + (d & 1);

    *sign = s & 1;
    *digit = d;
}

groupElement groupElement::operator*(const fieldElement &scalar) const {
    int skip;
    p224_felem nq[3], tmp[4];
    u64 bits;
    u32 sign, digit;

    p224_felem p_pre_comp[17][3];
    makePrecomp(p_pre_comp);
    /* set nq to the point at infinity */
    OPENSSL_memset(nq, 0, 3 * sizeof(p224_felem));

    /*
     * Loop over all scalars msb-to-lsb, interleaving additions of multiples
     * of the generator (two in each of the last 28 rounds) and additions of
     * other points multiples (every 5th round).
     */
    skip = 1;                   /* save two point operations in the first
                                 * round */
    size_t i;
    for (i = 220; i < 221; --i) {
        /* double */
        if (!skip)
            p224_point_double(nq[0], nq[1], nq[2], nq[0], nq[1], nq[2]);

        /* do other additions every 5 doublings */
        p224_felem_bytearray nbytes;
        scalar.toBin28(nbytes);
        if ((i % 5 == 0)) {
            /* loop over all scalars */
            bits = scalar.getBit(i + 4) << 5;
            bits |= scalar.getBit(i + 3) << 4;
            bits |= scalar.getBit(i + 2) << 3;
            bits |= scalar.getBit(i + 1) << 2;
            bits |= scalar.getBit(i) << 1;
            bits |= scalar.getBit(i - 1);

            ec_GFp_nistp_recode_scalar_bits(&sign, &digit, bits);

            /* select the point to add or subtract */
            p224_select_point(digit, 17, (const p224_felem(*)[3])p_pre_comp, tmp);
            p224_felem_neg(tmp[3], tmp[1]); /* (X, -Y, Z) is the negative point */
            p224_copy_conditional(tmp[1], tmp[3], sign);

            if (!skip) {
                p224_point_add(nq[0], nq[1], nq[2],
                          nq[0], nq[1], nq[2],
                          0, tmp[0], tmp[1], tmp[2]);
            } else {
                OPENSSL_memcpy(nq, tmp, 3 * sizeof(p224_felem));
                skip = 0;
            }
        }
    }

    groupElement res;
    p224_felem_contract(res.x, nq[0]);
    p224_felem_contract(res.y, nq[1]);
    p224_felem_contract(res.z, nq[2]);
    return res;
}

bool groupElement::operator==(const groupElement &other) const {
    fieldElement tmp1, tmp2, za23, zb23;
    fieldElement xa(x), ya(y), za(z), xb(other.x), yb(other.y), zb(other.z);
    zb23 = zb.sqr();
    tmp1 = xa * zb23;
    za23 = za.sqr();
    tmp2 = xb * za23;

    tmp1 -= tmp2;
    const u64 x_not_equal = !tmp1.isZero();

    zb23 *= zb;
    tmp1 = ya * zb23;
    za23 *= za;
    tmp2 = yb * za23;
    tmp1 -= tmp2;
    const bool y_not_equal = !tmp1.isZero();
    const bool x_and_y_equal = !(x_not_equal || y_not_equal);

    const bool a_not_infinity = !za.isZero();
    const bool b_not_infinity = !zb.isZero();
    const bool a_and_b_infinity = !(a_not_infinity || b_not_infinity);

    const bool equal =
            a_and_b_infinity || (a_not_infinity && b_not_infinity & x_and_y_equal);
    return equal;
}

bool groupElement::operator!=(const groupElement &other) const {
    return !(*this == other);
}

groupElement &groupElement::operator*=(const fieldElement &other) {
    *this = *this * other;
    return *this;
}

groupElement groupElement::dbl() const {
    groupElement res;
    p224_point_double(res.x, res.y, res.z, x, y, z);
    p224_felem_contract(res.x, res.x);
    p224_felem_contract(res.y, res.y);
    p224_felem_contract(res.z, res.z);
    return res;
}

groupElement groupElement::inv() const {
    groupElement res = *this;
    res.setInv();
    return res;
}

void groupElement::setDbl() {
    p224_point_double(x, y, z, x, y, z);
    p224_felem_contract(x, x);
    p224_felem_contract(y, y);
    p224_felem_contract(z, z);
}

void groupElement::setInv() {
    p224_felem y_prime;
    p224_felem_neg(y_prime, y);
    p224_felem_assign(y, y_prime);
}

void groupElement::makePrecomp(p224_felem (*out)[3]) const {
    OPENSSL_memset(out[0], 0, sizeof(p224_felem) * 3);

    p224_felem_assign(out[1][0], x);
    p224_felem_assign(out[1][1], y);
    p224_felem_assign(out[1][2], z);

    for (size_t j = 2; j <= 16; ++j) {
        if (j & 1) {
            p224_point_add(out[j][0], out[j][1], out[j][2], out[1][0], out[1][1],
                           out[1][2], 0, out[j - 1][0], out[j - 1][1], out[j - 1][2]);
        } else {
            p224_point_double(out[j][0], out[j][1], out[j][2], out[j >> 1][0],
                              out[j >> 1][1], out[j >> 1][2]);
        }
    }
}

groupElement::groupElement(const p224_felem _x, const p224_felem _y, const p224_felem _z) {
    OPENSSL_memcpy(x, _x, sizeof (p224_felem));
    OPENSSL_memcpy(y, _y, sizeof (p224_felem));
    OPENSSL_memcpy(z, _z, sizeof (p224_felem));
}

groupElement::groupElement(const p224_felem _p[3]) {
    OPENSSL_memcpy(x, _p[0], sizeof(p224_felem));
    OPENSSL_memcpy(y, _p[1], sizeof(p224_felem));
    OPENSSL_memcpy(z, _p[2], sizeof(p224_felem));
}

void groupElement::setInfinity() {
    OPENSSL_memset(x, 0, sizeof(p224_felem));
    OPENSSL_memset(y, 0, sizeof(p224_felem));
    OPENSSL_memset(z, 0, sizeof(p224_felem));
}

void groupElement::print() const {
    fprintf(stderr, "x: ");
    fieldElement(x).print(stderr);
    fprintf(stderr, "y: ");
    fieldElement(y).print(stderr);
    fprintf(stderr, "z: ");
    fieldElement(z).print(stderr);
}

bool groupElement::isOnCurve() {
    fieldElement xx(x), yy(y), zz(z), b(__B);

    // rh := X^2
    auto rh = xx.sqr();

    auto tmp = zz.sqr();
    auto z4 = tmp.sqr();
    auto z6 = z4 * tmp;

    // rh := rh + a*Z^4, a = 3
    tmp = z4 + z4;
    tmp = tmp + z4;
    rh = rh - tmp;

    // rh := (rh + a*Z^4)*X
    rh = rh * xx;
    // rh := rh + b*Z^6
    tmp = b * z6;
    rh = rh + tmp;

    // 'lh' := Y^2
    tmp = yy.sqr();

    tmp = tmp - rh;

    u64 mask = 0;
    for (int i = 0; i < 4; i++) {
        mask |= tmp[i];
    }

    auto not_equal = (bool) mask;

    // If Z = 0, the point is infinity, which is always on the curve.
    mask = 0;
    for (int i = 0; i < 4; i++) {
        mask |= zz[i];
    }

    auto not_infinity = (bool) mask;
    return !(not_infinity && not_equal);
}

}