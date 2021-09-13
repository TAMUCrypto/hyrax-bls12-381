//
// Created by juzix on 2021/6/9.
//

#ifndef HYRAX_P224_P224_HPP
#define HYRAX_P224_P224_HPP

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

#include <cstring>
#include "typedef.hpp"

extern const p224_felem __gens[3];
extern const p224_felem __A;
extern const p224_felem __B;

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

void *OPENSSL_memcpy(void *dst, const void *src, size_t n);

void *OPENSSL_memset(void *dst, int c, size_t n);
//
//void *OPENSSL_malloc(u64 sz);
//
//void OPENSSL_free(void *orig_ptr);

// Field element represented as a byte arrary. 28*8 = 224 bits is also the
// group order size for the elliptic curve, and we also use this type for
// scalars for point multiplication.
typedef u8 p224_felem_bytearray[28];

void p224_felem_one(p224_felem out);

unsigned char myrandom(void);

u64 p224_load_u64(const u8 in[8]);

// Helper functions to convert field elements to/from internal representation
void p224_bin28_to_felem(p224_felem out, const u8 in[28]);

void p224_felem_to_bin28(u8 out[28], const p224_felem in);


// Field operations, using the internal representation of field elements.
// NB! These operations are specific to our point multiplication and cannot be
// expected to be correct in general - e.g., multiplication with a large scalar
// will cause an overflow.

void p224_felem_assign(p224_felem out, const p224_felem in);

// Sum two field elements: out += in
void p224_felem_sum(p224_felem out, const p224_felem in);

// Sum with carry
void p224_felem_mysum(p224_felem out, const p224_felem in);

// Subtract field elements: out -= in
// Assumes in[i] < 2^57
void p224_felem_diff(p224_felem out, const p224_felem in);

// Diff with carry
void p224_felem_mydiff(p224_felem out, const p224_felem in);

// Subtract in unreduced 128-bit mode: out -= in
// Assumes in[i] < 2^119
void p224_widefelem_diff(p224_widefelem out, const p224_widefelem in);

// Subtract in mixed mode: out128 -= in64
// in[i] < 2^63
void p224_felem_diff_128_64(p224_widefelem out, const p224_felem in);

// Multiply a field element by a scalar: out = out * scalar
// The scalars we actually use are small, so results fit without overflow
void p224_felem_scalar(p224_felem out, const p224_limb scalar);

// Multiply an unreduced field element by a scalar: out = out * scalar
// The scalars we actually use are small, so results fit without overflow
void p224_widefelem_scalar(p224_widefelem out,
                                  const p224_widelimb scalar);

// Square a field element: out = in^2
void p224_felem_square(p224_widefelem out, const p224_felem in);

// Multiply two field elements: out = in1 * in2
void p224_felem_mul(p224_widefelem out, const p224_felem in1,
                           const p224_felem in2);

// Reduce seven 128-bit coefficients to four 64-bit coefficients.
// Requires in[i] < 2^126,
// ensures out[0] < 2^56, out[1] < 2^56, out[2] < 2^56, out[3] <= 2^56 + 2^16
void p224_felem_reduce(p224_felem out, const p224_widefelem in);

// Get negative value: out = -in
// Requires in[i] < 2^63,
// ensures out[0] < 2^56, out[1] < 2^56, out[2] < 2^56, out[3] <= 2^56 + 2^16
void p224_felem_neg(p224_felem out, const p224_felem in);

// Zero-check: returns 1 if input is 0, and 0 otherwise. We know that field
// elements are reduced to in < 2^225, so we only need to check three cases: 0,
// 2^224 - 2^96 + 1, and 2^225 - 2^97 + 2
p224_limb p224_felem_is_zero(const p224_felem in);

// Invert a field element
// Computation chain copied from djb's code
void p224_felem_inv(p224_felem out, const p224_felem in);

// Copy in constant time:
// if icopy == 1, copy in to out,
// if icopy == 0, copy out to itself.
void p224_copy_conditional(p224_felem out, const p224_felem in,
                                  p224_limb icopy);

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
                              const p224_felem y_in, const p224_felem z_in);

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
                           const p224_felem z2);

// p224_select_point selects the |idx|th point from a precomputation table and
// copies it to out.
void p224_select_point(const u64 idx, size_t size,
                              const p224_felem pre_comp[/*size*/][3],
                              p224_felem out[3]);

// p224_get_bit_from_felem returns the |i|th bit in |in|
u64 p224_get_bit(const p224_felem_bytearray in, size_t i);

u8 p224_get_bit_from_felem(const p224_limb *in, size_t bits);

// Reduce to unique minimal representation. Requires 0 <= in < 2*p (always
// call felem_reduce first)
void p224_felem_contract(p224_felem out, const p224_felem in);

bool p224_felem_equal(const p224_felem a, const p224_felem b);

bool p224_felem_non_zero_mask(const p224_limb *a);



#endif //HYRAX_P224_P224_HPP