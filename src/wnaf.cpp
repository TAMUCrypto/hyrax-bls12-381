//
// Created by juzix on 2021/6/11.
//

#include "wnaf.hpp"
#include "p224.hpp"

/* Originally written by Bodo Moeller for the OpenSSL project.
 * ====================================================================
 * Copyright (c) 1998-2005 The OpenSSL Project.  All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in
 *    the documentation and/or other materials provided with the
 *    distribution.
 *
 * 3. All advertising materials mentioning features or use of this
 *    software must display the following acknowledgment:
 *    "This product includes software developed by the OpenSSL Project
 *    for use in the OpenSSL Toolkit. (http://www.openssl.org/)"
 *
 * 4. The names "OpenSSL Toolkit" and "OpenSSL Project" must not be used to
 *    endorse or promote products derived from this software without
 *    prior written permission. For written permission, please contact
 *    openssl-core@openssl.org.
 *
 * 5. Products derived from this software may not be called "OpenSSL"
 *    nor may "OpenSSL" appear in their names without prior written
 *    permission of the OpenSSL Project.
 *
 * 6. Redistributions of any form whatsoever must retain the following
 *    acknowledgment:
 *    "This product includes software developed by the OpenSSL Project
 *    for use in the OpenSSL Toolkit (http://www.openssl.org/)"
 *
 * THIS SOFTWARE IS PROVIDED BY THE OpenSSL PROJECT ``AS IS'' AND ANY
 * EXPRESSED OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR __A PARTICULAR
 * PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE OpenSSL PROJECT OR
 * ITS CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT
 * NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
 * HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
 * STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED
 * OF THE POSSIBILITY OF SUCH DAMAGE.
 * ====================================================================
 *
 * This product includes cryptographic software written by Eric Young
 * (eay@cryptsoft.com).  This product includes software written by Tim
 * Hudson (tjh@cryptsoft.com).
 *
 */
/* ====================================================================
 * Copyright 2002 Sun Microsystems, Inc. ALL RIGHTS RESERVED.
 *
 * Portions of the attached software ("Contribution") are developed by
 * SUN MICROSYSTEMS, INC., and are contributed to the OpenSSL project.
 *
 * The Contribution is licensed pursuant to the OpenSSL open source
 * license provided above.
 *
 * The elliptic curve binary polynomial software is originally written by
 * Sheueling Chang Shantz and Douglas Stebila of Sun Microsystems
 * Laboratories. */

#include <assert.h>
#include <string.h>
#include "typedef.hpp"
#include "err.hpp"

// This file implements the wNAF-based interleaving multi-exponentiation method
// at:
//   http://link.springer.com/chapter/10.1007%2F3-540-45537-X_13
//   http://www.bmoeller.de/pdf/TI-01-08.multiexp.pdf

namespace hyrax_p224 {

namespace wnaf {
    vector<vector<p224_felem>> g_p224_pre_comp;

    void ec_compute_wNAF(i8 *out, const fieldElement &scalar, size_t bits, int w) {
        // 'i8' can represent integers with absolute values less than 2^7.
//    printf("bits: %d\n", bits);
        assert(0 < w && w <= 7);
        assert(bits != 0);
        int bit = 1 << w;         // 2^w, at most 128
        int next_bit = bit << 1;  // 2^(w+1), at most 256
        int mask = next_bit - 1;  // at most 255

        int window_val = scalar[0] & mask;
        for (size_t j = 0; j < bits + 1; j++) {
            assert(0 <= window_val && window_val <= next_bit);
            int digit = 0;
            if (window_val & 1) {
                assert(0 < window_val && window_val < next_bit);
                if (window_val & bit) {
                    digit = window_val - next_bit;
                    // We know -next_bit < digit < 0 and window_val - digit = next_bit.

                    // modified wNAF
                    if (j + w + 1 >= bits) {
                        // special case for generating modified wNAFs:
                        // no new bits will be added into window_val,
                        // so using a positive digit here will decrease
                        // the total length of the representation

                        digit = window_val & (mask >> 1);
                        // We know 0 < digit < bit and window_val - digit = bit.
                    }
                } else {
                    digit = window_val;
                    // We know 0 < digit < bit and window_val - digit = 0.
                }

                window_val -= digit;

                // Now window_val is 0 or 2^(w+1) in standard wNAF generation.
                // For modified window NAFs, it may also be 2^w.
                //
                // See the comments above for the derivation of each of these bounds.
                assert(window_val == 0 || window_val == next_bit || window_val == bit);
                assert(-bit < digit && digit < bit);

                // window_val was odd, so digit is also odd.
                assert(digit & 1);
            }

            out[j] = digit;

            // Incorporate the next bit. Previously, |window_val| <= |next_bit|, so if
            // we shift and add at most one copy of |bit|, this will continue to hold
            // afterwards.
            window_val >>= 1;
            window_val +=
                    bit * scalar.getBit(j + w + 1);
            assert(window_val <= next_bit);
        }

        // bits + 1 entries should be sufficient to consume all bits.
        assert(window_val == 0);
    }

// compute_precomp sets |out[i]| to (2*i+1)*p, for i from 0 to |len|.
    void compute_precomp(groupElement *out, const groupElement p, size_t len) {
        out[0] = p;
        groupElement two_p = p.dbl();
        for (size_t i = 1; i < len; i++)
            out[i] = out[i - 1] + two_p;

    }

    void lookup_precomp(groupElement &out, const groupElement *precomp, int digit) {
        out = digit < 0 ? precomp[-digit >> 1].inv() : precomp[digit >> 1];
    }

#define OPENSSL_ARRAY_SIZE(array) (sizeof(array) / sizeof((array)[0]))

#define EC_MAX_BYTES 66
#define EC_MAX_WORDS ((EC_MAX_BYTES + BN_BYTES - 1) / BN_BYTES)

// EC_WNAF_WINDOW_BITS is the window size to use for |ec_GFp_mont_mul_public|.
#define EC_WNAF_WINDOW_BITS 4

// EC_WNAF_TABLE_SIZE is the table size to use for |ec_GFp_mont_mul_public|.
#define EC_WNAF_TABLE_SIZE (1 << (EC_WNAF_WINDOW_BITS - 1))

// EC_WNAF_STACK is the number of points worth of data to stack-allocate and
// avoid a malloc.
#define EC_WNAF_STACK 3

    void mul_batch(groupElement &r, const fieldElement *g_scalar, const groupElement *points, const fieldElement *scalars, size_t num) {
        size_t wNAF_max_len = 224;
        size_t wNAF_len_g, *wNAF_lens;
        i8 wNAF_stack[EC_WNAF_STACK][EC_MAX_BYTES * 8 + 1];
        i8 (*wNAF_alloc)[EC_MAX_BYTES * 8 + 1] = NULL;
        i8 (*wNAF)[EC_MAX_BYTES * 8 + 1];
        groupElement precomp_stack[EC_WNAF_STACK][EC_WNAF_TABLE_SIZE];
        groupElement (*precomp_alloc)[EC_WNAF_TABLE_SIZE] = NULL;
        groupElement (*precomp)[EC_WNAF_TABLE_SIZE];
        if (num <= EC_WNAF_STACK) {
            wNAF = wNAF_stack;
            precomp = precomp_stack;
        } else {
            if (num >= ((size_t) -1) / sizeof(wNAF_alloc[0]) ||
                num >= ((size_t) -1) / sizeof(precomp_alloc[0])) {
                free(wNAF_alloc);
                free(precomp_alloc);
                free(wNAF_lens);
                exit(ERR_R_OVERFLOW);
            }
            wNAF_lens = (size_t *) malloc(num * sizeof(size_t));
            wNAF_alloc = (i8(*)[EC_MAX_BYTES * 8 + 1]) malloc(num * sizeof(wNAF_alloc[0]));
            precomp_alloc = (groupElement(*)[EC_WNAF_TABLE_SIZE]) malloc(num * sizeof(precomp_alloc[0]));
            if (wNAF_alloc == NULL || precomp_alloc == NULL) {
                free(wNAF_alloc);
                free(precomp_alloc);
                free(wNAF_lens);
                exit(ERR_R_MALLOC_FAILURE);
            }
            wNAF = wNAF_alloc;
            precomp = precomp_alloc;
        }

        i8 g_wNAF[EC_MAX_BYTES * 8 + 1];
        groupElement g_precomp[EC_WNAF_TABLE_SIZE];
        const groupElement g(__gens);
        wNAF_max_len = 0;
        if (g_scalar != NULL) {
            size_t bits = g_scalar->getBitWidth();
//        printf("bits: %d\n", bits);
            wNAF_len_g = bits + 1;
            assert(wNAF_len_g <= OPENSSL_ARRAY_SIZE(g_wNAF));
            wNAF_max_len = std::max(wNAF_max_len, wNAF_len_g);
            ec_compute_wNAF(g_wNAF, *g_scalar, bits, EC_WNAF_WINDOW_BITS);
            compute_precomp(g_precomp, g, EC_WNAF_TABLE_SIZE);
        }

        for (size_t i = 0; i < num; i++) {
            size_t bits = scalars[i].getBitWidth();
//        printf("bits: %d\n", bits);
            wNAF_lens[i] = bits + 1;
            wNAF_max_len = std::max(wNAF_max_len, wNAF_lens[i]);
            assert(wNAF_lens[i] <= OPENSSL_ARRAY_SIZE(wNAF[i]));
            ec_compute_wNAF(wNAF[i], scalars[i], bits, EC_WNAF_WINDOW_BITS);
            compute_precomp(precomp[i], points[i], EC_WNAF_TABLE_SIZE);
        }

        groupElement tmp;
        int r_is_at_infinity = 1;
//    printf("%lld\n", wNAF_max_len);
        for (size_t k = wNAF_max_len - 1; k < wNAF_max_len; k--) {
//        printf("k: %d\n", k);
            if (!r_is_at_infinity)
                r.setDbl();

            if (g_scalar != NULL && k < wNAF_len_g && g_wNAF[k] != 0) {
                lookup_precomp(tmp, g_precomp, g_wNAF[k]);
                if (r_is_at_infinity) {
                    r = tmp;
                    r_is_at_infinity = 0;
                } else {
                    r = r + tmp;
                }
            }

            for (size_t i = 0; i < num; i++) {
                if (k < wNAF_lens[i] && wNAF[i][k]) {
                    lookup_precomp(tmp, precomp[i], wNAF[i][k]);
                    if (r_is_at_infinity) {
                        r = tmp;
                        r_is_at_infinity = 0;
                    } else {
                        r = r + tmp;
                    }
                }
            }
        }

        if (r_is_at_infinity) r.setInfinity();
        free(wNAF_alloc);
        free(precomp_alloc);
        free(wNAF_lens);
    }

    groupElement vartimeMultiscalarMul(const vector<fieldElement> &exp, const vector<groupElement> &bas) {
        return vartimeMultiscalarMul(exp.begin(), bas.begin(), exp.size());
    }

    groupElement vartimeMultiscalarMul(const vector<fieldElement>::const_iterator &exp,
                                       const vector<groupElement>::const_iterator &bas, u64 n) {
        groupElement res;
        mul_batch(res, NULL, bas.base(), exp.base(), n);
        return res;
    }
}

}
