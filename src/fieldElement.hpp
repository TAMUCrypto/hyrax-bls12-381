//
// Created by juzix on 2021/6/1.
//

#ifndef HYRAX_P224_field_element_HPP
#define HYRAX_P224_field_element_HPP

#include <bits/stdc++.h>
#include "typedef.hpp"
#include "p224.hpp"

using std::vector;
using std::cerr;
using std::endl;
using std::ostream;
using std::istream;

namespace hyrax_p224 {
class fieldElement {
public:
    fieldElement(): data {0, 0, 0, 0} {}

    fieldElement(const fieldElement &b):
            data {b.data[0], b.data[1], b.data[2], b.data[3]}{}

    fieldElement(long long x):
            data {0, 0, 0, 0} {
        if (x >= 0) {
            data[0] = x & (1ULL << 56) - 1;
            data[1] = x >> 56;
        } else *this = zero() - fieldElement(-x);
    }

    explicit fieldElement(const u8 *in);

    explicit fieldElement(const p224_felem in);

    fieldElement operator+(const fieldElement &other) const;

    fieldElement operator-(const fieldElement &other) const;

    fieldElement operator - () const;

    fieldElement operator*(const fieldElement &other) const;

    bool operator==(const fieldElement &other) const;

    bool operator!=(const fieldElement &other) const;

    fieldElement &operator=(const fieldElement &other);

    fieldElement &operator+=(const fieldElement &other);

    fieldElement &operator-=(const fieldElement &other);

    fieldElement &operator*=(const fieldElement &other);

    friend ostream &operator << (ostream &out, const fieldElement &data);

    p224_limb operator [] (u8 i) const;
    p224_limb &operator [] (u8 i);

    explicit operator bool () const;

    [[nodiscard]] bool isNegative() const;

    [[nodiscard]] u8 getBitWidth() const;

    [[nodiscard]] u8 getBit(unsigned int i) const;

    char *toString() const;

    [[nodiscard]] __int128_t toint128() const;
    [[nodiscard]] bool lessThan(const fieldElement &b) const;

    [[nodiscard]] bool isZero();

    [[nodiscard]] fieldElement abs() const;
    [[nodiscard]] fieldElement sqr() const;
    [[nodiscard]] fieldElement inv() const;
    void setAbs();
    void setSqr();
    void setInv();

    void toBin28(p224_felem_bytearray bytes) const;

    void print(FILE *fileno) const;

    static void init();
    static fieldElement maxWithZero(const fieldElement &a, const fieldElement &b);
    static fieldElement maxUnsigned(const fieldElement &a, const fieldElement &b);
    static fieldElement get_root_of_unity(int order); //return a root of unity with order 2^[order]
    static fieldElement random();

    static fieldElement zero();
    static fieldElement one();
    static fieldElement *generateRandomness(unsigned int size);
    static fieldElement innerProd(vector<fieldElement>::iterator a, vector<fieldElement>::iterator b, u64 n);

    static bool initialized;
    static int multCounter, addCounter;
    static bool isCounting;
    static bool isSumchecking;

protected:
    p224_felem data;
};
}


#endif //HYRAX_P224_field_element_HPP
