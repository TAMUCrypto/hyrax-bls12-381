//
// Created by juzix on 2021/6/9.
//

/* Copyright (c) 2015, Google Inc.
 *
 * Permission to use, copy, modify, and/or distribute this software for any
 * purpose with or without fee is hereby granted, provided that the above
 * copyright notice and this permission notice appear in all copies.
 *
 * THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
 * WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
 * MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY
 * SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
 * WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN ACTION
 * OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF OR IN
 * CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE. */

// __A 64-bit implementation of the NIST P-224 elliptic curve point multiplication
//
// Inspired by Daniel J. Bernstein's public domain nistp224 implementation
// and Adam Langley's public domain 64-bit C implementation of curve25519.

#include <string.h>
#include <cassert>
#include <malloc.h>

#include "p224.hpp"

const p224_felem __gens[3] = {{0x3280d6115c1d21, 0xc1d356c2112234, 0x7f321390b94a03, 0xb70e0cbd6bb4bf},
                              {0xd5819985007e34, 0x75a05a07476444, 0xfb4c22dfe6cd43, 0xbd376388b5f723},
                              {1,                0,                0,                0}};
const p224_felem __A = {0xfffffffffffffe, 0xfffeffffffffff, 0xffffffffffffff, 0xffffffffffffff};
const p224_felem __B = {0xb39432355ffb4, 0xb0b7d7bfd8ba27, 0xabf54132565044, 0xb4050a850c04b3};

void *OPENSSL_memcpy(void *dst, const void *src, size_t n) {
    if (n == 0) {
        return dst;
    }

    return memcpy(dst, src, n);
}

void *OPENSSL_memset(void *dst, int c, size_t n) {
    if (n == 0) {
        return dst;
    }

    return memset(dst, c, n);
}

// https://bugs.llvm.org/show_bug.cgi?id=37598
#if defined(__ELF__) && defined(__GNUC__)
#define WEAK_SYMBOL_FUNC(rettype, name, args) \
  rettype name args __attribute__((weak));
#else
#define WEAK_SYMBOL_FUNC(rettype, name, args) static rettype(*name) args = NULL;
#endif

#if defined(OPENSSL_ASAN)
void __asan_poison_memory_region(const volatile void *addr, size_t size);
void __asan_unpoison_memory_region(const volatile void *addr, size_t size);
#else
void __asan_poison_memory_region(const void *addr, size_t size) {}
void __asan_unpoison_memory_region(const void *addr, size_t size) {}
#endif

// Field elements are represented as a_0 + 2^56*a_1 + 2^112*a_2 + 2^168*a_3
// using 64-bit coefficients called 'limbs', and sometimes (for multiplication
// results) as b_0 + 2^56*b_1 + 2^112*b_2 + 2^168*b_3 + 2^224*b_4 + 2^280*b_5 +
// 2^336*b_6 using 128-bit coefficients called 'widelimbs'. __A 4-p224_limb
// representation is an 'p224_felem'; a 7-p224_widelimb representation is a
// 'p224_widefelem'. Even within felems, bits of adjacent limbs overlap, and we
// don't always reduce the representations: we ensure that inputs to each
// p224_felem multiplication satisfy a_i < 2^60, so outputs satisfy b_i <
// 4*2^60*2^60, and fit into a 128-bit word without overflow. The coefficients
// are then again partially reduced to obtain an p224_felem satisfying a_i <
// 2^57. We only reduce to the unique minimal representation at the end of the
// computation.

// Precomputed multiples of the standard generator
// Points are given in coordinates (X, Y, Z) where Z normally is 1
// (0 for the point at infinity).
// For each field element, slice a_0 is word 0, etc.
//
// The table has 2 * 16 elements, starting with the following:
// index | bits    | point
// ------+---------+------------------------------
//     0 | 0 0 0 0 | 0G
//     1 | 0 0 0 1 | 1G
//     2 | 0 0 1 0 | 2^56G
//     3 | 0 0 1 1 | (2^56 + 1)__gens
//     4 | 0 1 0 0 | 2^112G
//     5 | 0 1 0 1 | (2^112 + 1)__gens
//     6 | 0 1 1 0 | (2^112 + 2^56)__gens
//     7 | 0 1 1 1 | (2^112 + 2^56 + 1)__gens
//     8 | 1 0 0 0 | 2^168G
//     9 | 1 0 0 1 | (2^168 + 1)__gens
//    10 | 1 0 1 0 | (2^168 + 2^56)__gens
//    11 | 1 0 1 1 | (2^168 + 2^56 + 1)__gens
//    12 | 1 1 0 0 | (2^168 + 2^112)__gens
//    13 | 1 1 0 1 | (2^168 + 2^112 + 1)__gens
//    14 | 1 1 1 0 | (2^168 + 2^112 + 2^56)__gens
//    15 | 1 1 1 1 | (2^168 + 2^112 + 2^56 + 1)__gens
// followed by a copy of this with each element multiplied by 2^28.
//
// The reason for this is so that we can clock bits into four different
// locations when doing simple scalar multiplies against the base point,
// and then another four locations using the second 16 elements.

static unsigned long in[4] = {0,0,0,0};
static unsigned long t[16];

#define ROTATE(x,b) x = (((x) << (b)) | ((x) >> (32 - (b))))

unsigned char myrandom(void) {
    static int pos = 0;
    unsigned long y;
    unsigned long x;
    unsigned long sum;
    int rounds;

    if (!pos) {
        x = 0;
        sum = 0;
        rounds = 3;
        t[1] =t[2] =t[3] =0; t[0] =in[0];
        t[5] =t[6] =t[7] =0; t[4] =in[1];
        t[9] =t[10]=t[11]=0; t[8] =in[2];
        t[13]=t[14]=t[15]=0; t[12]=in[3];
        do {
            sum += 0x9e3779b9; y=x; y+=sum;
            ROTATE(x,5);  x^=y; y=t[0];  x+=y; y=sum; t[0] =x; y+=x;
            ROTATE(x,7);  x^=y; y=t[1];  x+=y; y=sum; t[1] =x; y+=x;
            ROTATE(x,9);  x^=y; y=t[2];  x+=y; y=sum; t[2] =x; y+=x;
            ROTATE(x,13); x^=y; y=t[3];  x+=y; y=sum; t[3] =x; y+=x;
            ROTATE(x,5);  x^=y; y=t[4];  x+=y; y=sum; t[4] =x; y+=x;
            ROTATE(x,7);  x^=y; y=t[5];  x+=y; y=sum; t[5] =x; y+=x;
            ROTATE(x,9);  x^=y; y=t[6];  x+=y; y=sum; t[6] =x; y+=x;
            ROTATE(x,13); x^=y; y=t[7];  x+=y; y=sum; t[7] =x; y+=x;
            ROTATE(x,5);  x^=y; y=t[8];  x+=y; y=sum; t[8] =x; y+=x;
            ROTATE(x,7);  x^=y; y=t[9];  x+=y; y=sum; t[9] =x; y+=x;
            ROTATE(x,9);  x^=y; y=t[10]; x+=y; y=sum; t[10]=x; y+=x;
            ROTATE(x,13); x^=y; y=t[11]; x+=y; y=sum; t[11]=x; y+=x;
            ROTATE(x,5);  x^=y; y=t[12]; x+=y; y=sum; t[12]=x; y+=x;
            ROTATE(x,7);  x^=y; y=t[13]; x+=y; y=sum; t[13]=x; y+=x;
            ROTATE(x,9);  x^=y; y=t[14]; x+=y; y=sum; t[14]=x; y+=x;
            ROTATE(x,13); x^=y; y=t[15]; x+=y;
            t[15] = x;
        } while (--rounds);
        if (!++in[0]) if (!++in[1]) if (!++in[2]) ++in[3];
        pos = 16;
    }

    return t[--pos];
}

u64 p224_load_u64(const u8 in[8]) {
    u64 ret;
    OPENSSL_memcpy(&ret, in, sizeof(ret));
    return ret;
}

// Helper functions to convert field elements to/from internal representation
void p224_bin28_to_felem(p224_felem out, const u8 in[28]) {
    out[0] = p224_load_u64(in) & 0x00ffffffffffffff;
    out[1] = p224_load_u64(in + 7) & 0x00ffffffffffffff;
    out[2] = p224_load_u64(in + 14) & 0x00ffffffffffffff;
    out[3] = p224_load_u64(in + 20) >> 8;
}

void p224_felem_to_bin28(u8 out[28], const p224_felem in) {
    for (size_t i = 0; i < 7; ++i) {
        out[i] = in[0] >> (8 * i);
        out[i + 7] = in[1] >> (8 * i);
        out[i + 14] = in[2] >> (8 * i);
        out[i + 21] = in[3] >> (8 * i);
    }
}

// Requires 0 <= in < 2*p (always call p224_felem_reduce first)
void p224_felem_contract(p224_felem out, const p224_felem in) {
    // Reduce to unique minimal representation.
    static const i64 two56 = ((p224_limb)1) << 56;
    // 0 <= in < 2*p, p = 2^224 - 2^96 + 1
    // if in > p , reduce in = in - 2^224 + 2^96 - 1
    i64 tmp[4], a;
    tmp[0] = in[0];
    tmp[1] = in[1];
    tmp[2] = in[2];
    tmp[3] = in[3];
    // Case 1: a = 1 iff in >= 2^224
    a = (in[3] >> 56);
    tmp[0] -= a;
    tmp[1] += a << 40;
    tmp[3] &= 0x00ffffffffffffff;
    // Case 2: a = 0 iff p <= in < 2^224, i.e., the high 128 bits are all 1 and
    // the lower part is non-zero
    a = ((in[3] & in[2] & (in[1] | 0x000000ffffffffff)) + 1) |
        (((i64)(in[0] + (in[1] & 0x000000ffffffffff)) - 1) >> 63);
    a &= 0x00ffffffffffffff;
    // turn a into an all-one mask (if a = 0) or an all-zero mask
    a = (a - 1) >> 63;
    // subtract 2^224 - 2^96 + 1 if a is all-one
    tmp[3] &= a ^ 0xffffffffffffffff;
    tmp[2] &= a ^ 0xffffffffffffffff;
    tmp[1] &= (a ^ 0xffffffffffffffff) | 0x000000ffffffffff;
    tmp[0] -= 1 & a;

    // eliminate negative coefficients: if tmp[0] is negative, tmp[1] must
    // be non-zero, so we only need one step
    a = tmp[0] >> 63;
    tmp[0] += two56 & a;
    tmp[1] -= 1 & a;

    // carry 1 -> 2 -> 3
    tmp[2] += tmp[1] >> 56;
    tmp[1] &= 0x00ffffffffffffff;

    tmp[3] += tmp[2] >> 56;
    tmp[2] &= 0x00ffffffffffffff;

    // Now 0 <= tmp < p
    out[0] = tmp[0];
    out[1] = tmp[1];
    out[2] = tmp[2];
    out[3] = tmp[3];
}

void p224_felem_one(p224_felem out) {
    out[0] = 1;
    out[1] = 0;
    out[2] = 0;
    out[3] = 0;
}

// Field operations, using the internal representation of field elements.
// NB! These operations are specific to our point multiplication and cannot be
// expected to be correct in general - e.g., multiplication with a large scalar
// will cause an overflow.

void p224_felem_assign(p224_felem out, const p224_felem in) {
    out[0] = in[0];
    out[1] = in[1];
    out[2] = in[2];
    out[3] = in[3];
}

// Sum two field elements: out += in
void p224_felem_sum(p224_felem out, const p224_felem in) {
    out[0] += in[0];
    out[1] += in[1];
    out[2] += in[2];
    out[3] += in[3];
}

void p224_felem_mysum(p224_felem out, const p224_felem in) {
    out[0] += in[0];
    out[1] += in[1] + (out[0] >> 56);
    out[2] += in[2] + (out[1] >> 56);
    out[3] += in[3] + (out[2] >> 56);

    out[0] &= 0x00ffffffffffffff;
    out[1] &= 0x00ffffffffffffff;
    out[2] &= 0x00ffffffffffffff;
}

// Subtract field elements: out -= in
// Assumes in[i] < 2^57
void p224_felem_diff(p224_felem out, const p224_felem in) {
    static const p224_limb two58p2 =
            (((p224_limb)1) << 58) + (((p224_limb)1) << 2);
    static const p224_limb two58m2 =
            (((p224_limb)1) << 58) - (((p224_limb)1) << 2);
    static const p224_limb two58m42m2 =
            (((p224_limb)1) << 58) - (((p224_limb)1) << 42) - (((p224_limb)1) << 2);

    // Add 0 mod 2^224-2^96+1 to ensure out > in
    out[0] += two58p2;
    out[1] += two58m42m2;
    out[2] += two58m2;
    out[3] += two58m2;

    out[0] -= in[0];
    out[1] -= in[1];
    out[2] -= in[2];
    out[3] -= in[3];
}

void p224_felem_mydiff(p224_felem out, const p224_felem in) {
    static const p224_limb two58p2 =
            (((p224_limb)1) << 58) + (((p224_limb)1) << 2);
    static const p224_limb two58m2 =
            (((p224_limb)1) << 58) - (((p224_limb)1) << 2);
    static const p224_limb two58m42m2 =
            (((p224_limb)1) << 58) - (((p224_limb)1) << 42) - (((p224_limb)1) << 2);

    // Add 0 mod 2^224-2^96+1 to ensure out > in
    out[0] += two58p2;
    out[1] += two58m42m2;
    out[2] += two58m2;
    out[3] += two58m2;

    out[0] = out[0] - in[0] + (1ULL << 56);
    out[1] = out[1] - in[1] + 0x00ffffffffffffff + (out[0] >> 56);
    out[2] = out[2] - in[2] + 0x00ffffffffffffff + (out[1] >> 56);
    out[3] = out[3] - in[3] - 1 + (out[2] >> 56);

    out[0] &= 0x00ffffffffffffff;
    out[1] &= 0x00ffffffffffffff;
    out[2] &= 0x00ffffffffffffff;
}

// Subtract in unreduced 128-bit mode: out -= in
// Assumes in[i] < 2^119
void p224_widefelem_diff(p224_widefelem out, const p224_widefelem in) {
    static const p224_widelimb two120 = ((p224_widelimb)1) << 120;
    static const p224_widelimb two120m64 =
            (((p224_widelimb)1) << 120) - (((p224_widelimb)1) << 64);
    static const p224_widelimb two120m104m64 = (((p224_widelimb)1) << 120) -
                                               (((p224_widelimb)1) << 104) -
                                               (((p224_widelimb)1) << 64);

    // Add 0 mod 2^224-2^96+1 to ensure out > in
    out[0] += two120;
    out[1] += two120m64;
    out[2] += two120m64;
    out[3] += two120;
    out[4] += two120m104m64;
    out[5] += two120m64;
    out[6] += two120m64;

    out[0] -= in[0];
    out[1] -= in[1];
    out[2] -= in[2];
    out[3] -= in[3];
    out[4] -= in[4];
    out[5] -= in[5];
    out[6] -= in[6];
}

// Subtract in mixed mode: out128 -= in64
// in[i] < 2^63
void p224_felem_diff_128_64(p224_widefelem out, const p224_felem in) {
    static const p224_widelimb two64p8 =
            (((p224_widelimb)1) << 64) + (((p224_widelimb)1) << 8);
    static const p224_widelimb two64m8 =
            (((p224_widelimb)1) << 64) - (((p224_widelimb)1) << 8);
    static const p224_widelimb two64m48m8 = (((p224_widelimb)1) << 64) -
                                            (((p224_widelimb)1) << 48) -
                                            (((p224_widelimb)1) << 8);

    // Add 0 mod 2^224-2^96+1 to ensure out > in
    out[0] += two64p8;
    out[1] += two64m48m8;
    out[2] += two64m8;
    out[3] += two64m8;

    out[0] -= in[0];
    out[1] -= in[1];
    out[2] -= in[2];
    out[3] -= in[3];
}

// Multiply a field element by a scalar: out = out * scalar
// The scalars we actually use are small, so results fit without overflow
void p224_felem_scalar(p224_felem out, const p224_limb scalar) {
    out[0] *= scalar;
    out[1] *= scalar;
    out[2] *= scalar;
    out[3] *= scalar;
}

// Multiply an unreduced field element by a scalar: out = out * scalar
// The scalars we actually use are small, so results fit without overflow
void p224_widefelem_scalar(p224_widefelem out,
                                  const p224_widelimb scalar) {
    out[0] *= scalar;
    out[1] *= scalar;
    out[2] *= scalar;
    out[3] *= scalar;
    out[4] *= scalar;
    out[5] *= scalar;
    out[6] *= scalar;
}

// Square a field element: out = in^2
void p224_felem_square(p224_widefelem out, const p224_felem in) {
    p224_limb tmp0, tmp1, tmp2;
    tmp0 = 2 * in[0];
    tmp1 = 2 * in[1];
    tmp2 = 2 * in[2];
    out[0] = ((p224_widelimb)in[0]) * in[0];
    out[1] = ((p224_widelimb)in[0]) * tmp1;
    out[2] = ((p224_widelimb)in[0]) * tmp2 + ((p224_widelimb)in[1]) * in[1];
    out[3] = ((p224_widelimb)in[3]) * tmp0 + ((p224_widelimb)in[1]) * tmp2;
    out[4] = ((p224_widelimb)in[3]) * tmp1 + ((p224_widelimb)in[2]) * in[2];
    out[5] = ((p224_widelimb)in[3]) * tmp2;
    out[6] = ((p224_widelimb)in[3]) * in[3];
}

// Multiply two field elements: out = in1 * in2
void p224_felem_mul(p224_widefelem out, const p224_felem in1,
                           const p224_felem in2) {
    out[0] = ((p224_widelimb)in1[0]) * in2[0];
    out[1] = ((p224_widelimb)in1[0]) * in2[1] + ((p224_widelimb)in1[1]) * in2[0];
    out[2] = ((p224_widelimb)in1[0]) * in2[2] + ((p224_widelimb)in1[1]) * in2[1] +
             ((p224_widelimb)in1[2]) * in2[0];
    out[3] = ((p224_widelimb)in1[0]) * in2[3] + ((p224_widelimb)in1[1]) * in2[2] +
             ((p224_widelimb)in1[2]) * in2[1] + ((p224_widelimb)in1[3]) * in2[0];
    out[4] = ((p224_widelimb)in1[1]) * in2[3] + ((p224_widelimb)in1[2]) * in2[2] +
             ((p224_widelimb)in1[3]) * in2[1];
    out[5] = ((p224_widelimb)in1[2]) * in2[3] + ((p224_widelimb)in1[3]) * in2[2];
    out[6] = ((p224_widelimb)in1[3]) * in2[3];
}

// Reduce seven 128-bit coefficients to four 64-bit coefficients.
// Requires in[i] < 2^126,
// ensures out[0] < 2^56, out[1] < 2^56, out[2] < 2^56, out[3] <= 2^56 + 2^16
void p224_felem_reduce(p224_felem out, const p224_widefelem in) {
    static const p224_widelimb two127p15 =
            (((p224_widelimb)1) << 127) + (((p224_widelimb)1) << 15);
    static const p224_widelimb two127m71 =
            (((p224_widelimb)1) << 127) - (((p224_widelimb)1) << 71);
    static const p224_widelimb two127m71m55 = (((p224_widelimb)1) << 127) -
                                              (((p224_widelimb)1) << 71) -
                                              (((p224_widelimb)1) << 55);
    p224_widelimb output[5];

    // Add 0 mod 2^224-2^96+1 to ensure all differences are positive
    output[0] = in[0] + two127p15;
    output[1] = in[1] + two127m71m55;
    output[2] = in[2] + two127m71;
    output[3] = in[3];
    output[4] = in[4];

    // Eliminate in[4], in[5], in[6]
    output[4] += in[6] >> 16;
    output[3] += (in[6] & 0xffff) << 40;
    output[2] -= in[6];

    output[3] += in[5] >> 16;
    output[2] += (in[5] & 0xffff) << 40;
    output[1] -= in[5];

    output[2] += output[4] >> 16;
    output[1] += (output[4] & 0xffff) << 40;
    output[0] -= output[4];

    // Carry 2 -> 3 -> 4
    output[3] += output[2] >> 56;
    output[2] &= 0x00ffffffffffffff;

    output[4] = output[3] >> 56;
    output[3] &= 0x00ffffffffffffff;

    // Now output[2] < 2^56, output[3] < 2^56, output[4] < 2^72

    // Eliminate output[4]
    output[2] += output[4] >> 16;
    // output[2] < 2^56 + 2^56 = 2^57
    output[1] += (output[4] & 0xffff) << 40;
    output[0] -= output[4];

    // Carry 0 -> 1 -> 2 -> 3
    output[1] += output[0] >> 56;
    out[0] = output[0] & 0x00ffffffffffffff;

    output[2] += output[1] >> 56;
    // output[2] < 2^57 + 2^72
    out[1] = output[1] & 0x00ffffffffffffff;
    output[3] += output[2] >> 56;
    // output[3] <= 2^56 + 2^16
    out[2] = output[2] & 0x00ffffffffffffff;

    // out[0] < 2^56, out[1] < 2^56, out[2] < 2^56,
    // out[3] <= 2^56 + 2^16 (due to final carry),
    // so out < 2*p
    out[3] = output[3];
}

// Get negative value: out = -in
// Requires in[i] < 2^63,
// ensures out[0] < 2^56, out[1] < 2^56, out[2] < 2^56, out[3] <= 2^56 + 2^16
void p224_felem_neg(p224_felem out, const p224_felem in) {
    p224_widefelem tmp = {0};
    p224_felem_diff_128_64(tmp, in);
    p224_felem_reduce(out, tmp);
}

// Zero-check: returns 1 if input is 0, and 0 otherwise. We know that field
// elements are reduced to in < 2^225, so we only need to check three cases: 0,
// 2^224 - 2^96 + 1, and 2^225 - 2^97 + 2
p224_limb p224_felem_is_zero(const p224_felem in) {
    p224_limb zero = in[0] | in[1] | in[2] | in[3];
    zero = (((i64)(zero)-1) >> 63) & 1;

    p224_limb two224m96p1 = (in[0] ^ 1) | (in[1] ^ 0x00ffff0000000000) |
                            (in[2] ^ 0x00ffffffffffffff) |
                            (in[3] ^ 0x00ffffffffffffff);
    two224m96p1 = (((i64)(two224m96p1)-1) >> 63) & 1;
    p224_limb two225m97p2 = (in[0] ^ 2) | (in[1] ^ 0x00fffe0000000000) |
                            (in[2] ^ 0x00ffffffffffffff) |
                            (in[3] ^ 0x01ffffffffffffff);
    two225m97p2 = (((i64)(two225m97p2)-1) >> 63) & 1;
    return (zero | two224m96p1 | two225m97p2);
}

// Invert a field element
// Computation chain copied from djb's code
void p224_felem_inv(p224_felem out, const p224_felem in) {
    p224_felem ftmp, ftmp2, ftmp3, ftmp4;
    p224_widefelem tmp;

    p224_felem_square(tmp, in);
    p224_felem_reduce(ftmp, tmp);  // 2
    p224_felem_mul(tmp, in, ftmp);
    p224_felem_reduce(ftmp, tmp);  // 2^2 - 1
    p224_felem_square(tmp, ftmp);
    p224_felem_reduce(ftmp, tmp);  // 2^3 - 2
    p224_felem_mul(tmp, in, ftmp);
    p224_felem_reduce(ftmp, tmp);  // 2^3 - 1
    p224_felem_square(tmp, ftmp);
    p224_felem_reduce(ftmp2, tmp);  // 2^4 - 2
    p224_felem_square(tmp, ftmp2);
    p224_felem_reduce(ftmp2, tmp);  // 2^5 - 4
    p224_felem_square(tmp, ftmp2);
    p224_felem_reduce(ftmp2, tmp);  // 2^6 - 8
    p224_felem_mul(tmp, ftmp2, ftmp);
    p224_felem_reduce(ftmp, tmp);  // 2^6 - 1
    p224_felem_square(tmp, ftmp);
    p224_felem_reduce(ftmp2, tmp);  // 2^7 - 2
    for (size_t i = 0; i < 5; ++i) {  // 2^12 - 2^6
        p224_felem_square(tmp, ftmp2);
        p224_felem_reduce(ftmp2, tmp);
    }
    p224_felem_mul(tmp, ftmp2, ftmp);
    p224_felem_reduce(ftmp2, tmp);  // 2^12 - 1
    p224_felem_square(tmp, ftmp2);
    p224_felem_reduce(ftmp3, tmp);  // 2^13 - 2
    for (size_t i = 0; i < 11; ++i) {  // 2^24 - 2^12
        p224_felem_square(tmp, ftmp3);
        p224_felem_reduce(ftmp3, tmp);
    }
    p224_felem_mul(tmp, ftmp3, ftmp2);
    p224_felem_reduce(ftmp2, tmp);  // 2^24 - 1
    p224_felem_square(tmp, ftmp2);
    p224_felem_reduce(ftmp3, tmp);  // 2^25 - 2
    for (size_t i = 0; i < 23; ++i) {  // 2^48 - 2^24
        p224_felem_square(tmp, ftmp3);
        p224_felem_reduce(ftmp3, tmp);
    }
    p224_felem_mul(tmp, ftmp3, ftmp2);
    p224_felem_reduce(ftmp3, tmp);  // 2^48 - 1
    p224_felem_square(tmp, ftmp3);
    p224_felem_reduce(ftmp4, tmp);  // 2^49 - 2
    for (size_t i = 0; i < 47; ++i) {  // 2^96 - 2^48
        p224_felem_square(tmp, ftmp4);
        p224_felem_reduce(ftmp4, tmp);
    }
    p224_felem_mul(tmp, ftmp3, ftmp4);
    p224_felem_reduce(ftmp3, tmp);  // 2^96 - 1
    p224_felem_square(tmp, ftmp3);
    p224_felem_reduce(ftmp4, tmp);  // 2^97 - 2
    for (size_t i = 0; i < 23; ++i) {  // 2^120 - 2^24
        p224_felem_square(tmp, ftmp4);
        p224_felem_reduce(ftmp4, tmp);
    }
    p224_felem_mul(tmp, ftmp2, ftmp4);
    p224_felem_reduce(ftmp2, tmp);  // 2^120 - 1
    for (size_t i = 0; i < 6; ++i) {  // 2^126 - 2^6
        p224_felem_square(tmp, ftmp2);
        p224_felem_reduce(ftmp2, tmp);
    }
    p224_felem_mul(tmp, ftmp2, ftmp);
    p224_felem_reduce(ftmp, tmp);  // 2^126 - 1
    p224_felem_square(tmp, ftmp);
    p224_felem_reduce(ftmp, tmp);  // 2^127 - 2
    p224_felem_mul(tmp, ftmp, in);
    p224_felem_reduce(ftmp, tmp);  // 2^127 - 1
    for (size_t i = 0; i < 97; ++i) {  // 2^224 - 2^97
        p224_felem_square(tmp, ftmp);
        p224_felem_reduce(ftmp, tmp);
    }
    p224_felem_mul(tmp, ftmp, ftmp3);
    p224_felem_reduce(out, tmp);  // 2^224 - 2^96 - 1
}

// Copy in constant time:
// if icopy == 1, copy in to out,
// if icopy == 0, copy out to itself.
void p224_copy_conditional(p224_felem out, const p224_felem in,
                                  p224_limb icopy) {
    // icopy is a (64-bit) 0 or 1, so copy is either all-zero or all-one
    const p224_limb copy = -icopy;
    for (size_t i = 0; i < 4; ++i) {
        const p224_limb tmp = copy & (in[i] ^ out[i]);
        out[i] ^= tmp;
    }
}

// ELLIPTIC CURVE POINT OPERATIONS
//
// Points are represented in Jacobian projective coordinates:
// (X, Y, Z) corresponds to the affine point (X/Z^2, Y/Z^3),
// or to the point at infinity if Z == 0.

// Double an elliptic curve point:
// (X', Y', Z') = 2 * (X, Y, Z), where
// X' = (3 * (X - Z^2) * (X + Z^2))^2 - 8 * X * Y^2
// Y' = 3 * (X - Z^2) * (X + Z^2) * (4 * X * Y^2 - X') - 8 * Y^2
// Z' = (Y + Z)^2 - Y^2 - Z^2 = 2 * Y * Z
// Outputs can equal corresponding inputs, i.e., x_out == x_in is allowed,
// while x_out == y_in is not (maybe this works, but it's not tested).
void p224_point_double(p224_felem x_out, p224_felem y_out,
                              p224_felem z_out, const p224_felem x_in,
                              const p224_felem y_in, const p224_felem z_in) {
    p224_widefelem tmp, tmp2;
    p224_felem delta, gamma, beta, alpha, ftmp, ftmp2;

    p224_felem_assign(ftmp, x_in);
    p224_felem_assign(ftmp2, x_in);

    // delta = z^2
    p224_felem_square(tmp, z_in);
    p224_felem_reduce(delta, tmp);

    // gamma = y^2
    p224_felem_square(tmp, y_in);
    p224_felem_reduce(gamma, tmp);

    // beta = x*gamma
    p224_felem_mul(tmp, x_in, gamma);
    p224_felem_reduce(beta, tmp);

    // alpha = 3*(x-delta)*(x+delta)
    p224_felem_diff(ftmp, delta);
    // ftmp[i] < 2^57 + 2^58 + 2 < 2^59
    p224_felem_sum(ftmp2, delta);
    // ftmp2[i] < 2^57 + 2^57 = 2^58
    p224_felem_scalar(ftmp2, 3);
    // ftmp2[i] < 3 * 2^58 < 2^60
    p224_felem_mul(tmp, ftmp, ftmp2);
    // tmp[i] < 2^60 * 2^59 * 4 = 2^121
    p224_felem_reduce(alpha, tmp);

    // x' = alpha^2 - 8*beta
    p224_felem_square(tmp, alpha);
    // tmp[i] < 4 * 2^57 * 2^57 = 2^116
    p224_felem_assign(ftmp, beta);
    p224_felem_scalar(ftmp, 8);
    // ftmp[i] < 8 * 2^57 = 2^60
    p224_felem_diff_128_64(tmp, ftmp);
    // tmp[i] < 2^116 + 2^64 + 8 < 2^117
    p224_felem_reduce(x_out, tmp);

    // z' = (y + z)^2 - gamma - delta
    p224_felem_sum(delta, gamma);
    // delta[i] < 2^57 + 2^57 = 2^58
    p224_felem_assign(ftmp, y_in);
    p224_felem_sum(ftmp, z_in);
    // ftmp[i] < 2^57 + 2^57 = 2^58
    p224_felem_square(tmp, ftmp);
    // tmp[i] < 4 * 2^58 * 2^58 = 2^118
    p224_felem_diff_128_64(tmp, delta);
    // tmp[i] < 2^118 + 2^64 + 8 < 2^119
    p224_felem_reduce(z_out, tmp);

    // y' = alpha*(4*beta - x') - 8*gamma^2
    p224_felem_scalar(beta, 4);
    // beta[i] < 4 * 2^57 = 2^59
    p224_felem_diff(beta, x_out);
    // beta[i] < 2^59 + 2^58 + 2 < 2^60
    p224_felem_mul(tmp, alpha, beta);
    // tmp[i] < 4 * 2^57 * 2^60 = 2^119
    p224_felem_square(tmp2, gamma);
    // tmp2[i] < 4 * 2^57 * 2^57 = 2^116
    p224_widefelem_scalar(tmp2, 8);
    // tmp2[i] < 8 * 2^116 = 2^119
    p224_widefelem_diff(tmp, tmp2);
    // tmp[i] < 2^119 + 2^120 < 2^121
    p224_felem_reduce(y_out, tmp);
}

// Add two elliptic curve points:
// (X_1, Y_1, Z_1) + (X_2, Y_2, Z_2) = (X_3, Y_3, Z_3), where
// X_3 = (Z_1^3 * Y_2 - Z_2^3 * Y_1)^2 - (Z_1^2 * X_2 - Z_2^2 * X_1)^3 -
// 2 * Z_2^2 * X_1 * (Z_1^2 * X_2 - Z_2^2 * X_1)^2
// Y_3 = (Z_1^3 * Y_2 - Z_2^3 * Y_1) * (Z_2^2 * X_1 * (Z_1^2 * X_2 - Z_2^2 *
// X_1)^2 - X_3) -
//        Z_2^3 * Y_1 * (Z_1^2 * X_2 - Z_2^2 * X_1)^3
// Z_3 = (Z_1^2 * X_2 - Z_2^2 * X_1) * (Z_1 * Z_2)
//
// This runs faster if 'mixed' is set, which requires Z_2 = 1 or Z_2 = 0.

// This function is not entirely constant-time: it includes a branch for
// checking whether the two input points are equal, (while not equal to the
// point at infinity). This case never happens during single point
// multiplication, so there is no timing leak for ECDH or ECDSA signing.
void p224_point_add(p224_felem x3, p224_felem y3, p224_felem z3,
                           const p224_felem x1, const p224_felem y1,
                           const p224_felem z1, const int mixed,
                           const p224_felem x2, const p224_felem y2,
                           const p224_felem z2) {

    p224_felem ftmp, ftmp2, ftmp3, ftmp4, ftmp5, x_out, y_out, z_out;
    p224_widefelem tmp, tmp2;
    p224_limb z1_is_zero, z2_is_zero, x_equal, y_equal;

    if (!mixed) {
        // ftmp2 = z2^2
        p224_felem_square(tmp, z2);
        p224_felem_reduce(ftmp2, tmp);

        // ftmp4 = z2^3
        p224_felem_mul(tmp, ftmp2, z2);
        p224_felem_reduce(ftmp4, tmp);

        // ftmp4 = z2^3*y1
        p224_felem_mul(tmp2, ftmp4, y1);
        p224_felem_reduce(ftmp4, tmp2);

        // ftmp2 = z2^2*x1
        p224_felem_mul(tmp2, ftmp2, x1);
        p224_felem_reduce(ftmp2, tmp2);
    } else {
        // We'll assume z2 = 1 (special case z2 = 0 is handled later)

        // ftmp4 = z2^3*y1
        p224_felem_assign(ftmp4, y1);

        // ftmp2 = z2^2*x1
        p224_felem_assign(ftmp2, x1);
    }

    // ftmp = z1^2
    p224_felem_square(tmp, z1);
    p224_felem_reduce(ftmp, tmp);

    // ftmp3 = z1^3
    p224_felem_mul(tmp, ftmp, z1);
    p224_felem_reduce(ftmp3, tmp);

    // tmp = z1^3*y2
    p224_felem_mul(tmp, ftmp3, y2);
    // tmp[i] < 4 * 2^57 * 2^57 = 2^116

    // ftmp3 = z1^3*y2 - z2^3*y1
    p224_felem_diff_128_64(tmp, ftmp4);
    // tmp[i] < 2^116 + 2^64 + 8 < 2^117
    p224_felem_reduce(ftmp3, tmp);

    // tmp = z1^2*x2
    p224_felem_mul(tmp, ftmp, x2);
    // tmp[i] < 4 * 2^57 * 2^57 = 2^116

    // ftmp = z1^2*x2 - z2^2*x1
    p224_felem_diff_128_64(tmp, ftmp2);
    // tmp[i] < 2^116 + 2^64 + 8 < 2^117
    p224_felem_reduce(ftmp, tmp);

    // the formulae are incorrect if the points are equal
    // so we check for this and do doubling if this happens
    x_equal = p224_felem_is_zero(ftmp);
    y_equal = p224_felem_is_zero(ftmp3);
    z1_is_zero = p224_felem_is_zero(z1);
    z2_is_zero = p224_felem_is_zero(z2);
    // In affine coordinates, (X_1, Y_1) == (X_2, Y_2)
    p224_limb is_nontrivial_double =
            x_equal & y_equal & (1 - z1_is_zero) & (1 - z2_is_zero);
    if (is_nontrivial_double) {
        p224_point_double(x3, y3, z3, x1, y1, z1);
        return;
    }

    // ftmp5 = z1*z2
    if (!mixed) {
        p224_felem_mul(tmp, z1, z2);
        p224_felem_reduce(ftmp5, tmp);
    } else {
        // special case z2 = 0 is handled later
        p224_felem_assign(ftmp5, z1);
    }

    // z_out = (z1^2*x2 - z2^2*x1)*(z1*z2)
    p224_felem_mul(tmp, ftmp, ftmp5);
    p224_felem_reduce(z_out, tmp);

    // ftmp = (z1^2*x2 - z2^2*x1)^2
    p224_felem_assign(ftmp5, ftmp);
    p224_felem_square(tmp, ftmp);
    p224_felem_reduce(ftmp, tmp);

    // ftmp5 = (z1^2*x2 - z2^2*x1)^3
    p224_felem_mul(tmp, ftmp, ftmp5);
    p224_felem_reduce(ftmp5, tmp);

    // ftmp2 = z2^2*x1*(z1^2*x2 - z2^2*x1)^2
    p224_felem_mul(tmp, ftmp2, ftmp);
    p224_felem_reduce(ftmp2, tmp);

    // tmp = z2^3*y1*(z1^2*x2 - z2^2*x1)^3
    p224_felem_mul(tmp, ftmp4, ftmp5);
    // tmp[i] < 4 * 2^57 * 2^57 = 2^116

    // tmp2 = (z1^3*y2 - z2^3*y1)^2
    p224_felem_square(tmp2, ftmp3);
    // tmp2[i] < 4 * 2^57 * 2^57 < 2^116

    // tmp2 = (z1^3*y2 - z2^3*y1)^2 - (z1^2*x2 - z2^2*x1)^3
    p224_felem_diff_128_64(tmp2, ftmp5);
    // tmp2[i] < 2^116 + 2^64 + 8 < 2^117

    // ftmp5 = 2*z2^2*x1*(z1^2*x2 - z2^2*x1)^2
    p224_felem_assign(ftmp5, ftmp2);
    p224_felem_scalar(ftmp5, 2);
    // ftmp5[i] < 2 * 2^57 = 2^58

    /* x_out = (z1^3*y2 - z2^3*y1)^2 - (z1^2*x2 - z2^2*x1)^3 -
       2*z2^2*x1*(z1^2*x2 - z2^2*x1)^2 */
    p224_felem_diff_128_64(tmp2, ftmp5);
    // tmp2[i] < 2^117 + 2^64 + 8 < 2^118
    p224_felem_reduce(x_out, tmp2);

    // ftmp2 = z2^2*x1*(z1^2*x2 - z2^2*x1)^2 - x_out
    p224_felem_diff(ftmp2, x_out);
    // ftmp2[i] < 2^57 + 2^58 + 2 < 2^59

    // tmp2 = (z1^3*y2 - z2^3*y1)*(z2^2*x1*(z1^2*x2 - z2^2*x1)^2 - x_out)
    p224_felem_mul(tmp2, ftmp3, ftmp2);
    // tmp2[i] < 4 * 2^57 * 2^59 = 2^118

    /* y_out = (z1^3*y2 - z2^3*y1)*(z2^2*x1*(z1^2*x2 - z2^2*x1)^2 - x_out) -
       z2^3*y1*(z1^2*x2 - z2^2*x1)^3 */
    p224_widefelem_diff(tmp2, tmp);
    // tmp2[i] < 2^118 + 2^120 < 2^121
    p224_felem_reduce(y_out, tmp2);

    // the result (x_out, y_out, z_out) is incorrect if one of the inputs is
    // the point at infinity, so we need to check for this separately

    // if point 1 is at infinity, copy point 2 to output, and vice versa
    p224_copy_conditional(x_out, x2, z1_is_zero);
    p224_copy_conditional(x_out, x1, z2_is_zero);
    p224_copy_conditional(y_out, y2, z1_is_zero);
    p224_copy_conditional(y_out, y1, z2_is_zero);
    p224_copy_conditional(z_out, z2, z1_is_zero);
    p224_copy_conditional(z_out, z1, z2_is_zero);
    p224_felem_assign(x3, x_out);
    p224_felem_assign(y3, y_out);
    p224_felem_assign(z3, z_out);
}

// p224_select_point selects the |idx|th point from a precomputation table and
// copies it to out.
void p224_select_point(const u64 idx, size_t size,
                              const p224_felem pre_comp[/*size*/][3],
                              p224_felem out[3]) {
    p224_limb *outlimbs = &out[0][0];
    OPENSSL_memset(outlimbs, 0, 3 * sizeof(p224_felem));

    for (size_t i = 0; i < size; i++) {
        const p224_limb *inlimbs = &pre_comp[i][0][0];
        u64 mask = i ^ idx;
        mask |= mask >> 4;
        mask |= mask >> 2;
        mask |= mask >> 1;
        mask &= 1;
        mask--;
        for (size_t j = 0; j < 4 * 3; j++) {
            outlimbs[j] |= inlimbs[j] & mask;
        }
    }
}

// p224_get_bit returns the |i|th bit in |in|
u64 p224_get_bit(const p224_felem_bytearray in, size_t i) {
    if (i >= 224) {
        return 0;
    }
    return (in[i >> 3] >> (i & 7)) & 1;
}

u8 p224_get_bit_from_felem(const p224_felem in, size_t bits) {
    if (bits >= 224) return 0;
    size_t i = bits / 56;
    size_t j = bits % 56;
    return (in[i] >> j) & 1;
}

bool p224_felem_equal(const p224_limb *a, const p224_limb *b) {
    return memcmp(a, b, sizeof(p224_felem)) == 0;
}

bool p224_felem_non_zero_mask(const p224_limb *a) {
    u64 mask = 0;
    for (int i = 0; i < 4; i++) {
        mask |= a[i];
    }
    return (bool) mask;
}
