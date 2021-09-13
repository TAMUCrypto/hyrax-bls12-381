//
// Created by juzix on 2021/6/1.
//

#ifndef HYRAX_P224_GROUPELEMENT_HPP
#define HYRAX_P224_GROUPELEMENT_HPP

#include <bits/stdc++.h>
#include <openssl/obj_mac.h>
#include <openssl/ec.h>
#include <openssl/bn.h>
#include "typedef.hpp"
#include "fieldElement.hpp"

using std::vector;

namespace hyrax_p224 {

class groupElement {
public:
    groupElement(const p224_felem _x, const p224_felem _y, const p224_felem _z);

    explicit groupElement(const p224_felem _p[3]);

    groupElement() = default;

    static groupElement random();

    static groupElement zero();

    [[nodiscard]] groupElement dbl() const;
    [[nodiscard]] groupElement inv() const;

    void setDbl();
    void setInv();
    void setInfinity();

    groupElement operator + (const groupElement &other) const;
    groupElement operator * (const fieldElement &other) const;
    bool operator == (const groupElement &other) const;
    bool operator != (const groupElement &other) const;

    groupElement &operator = (const groupElement &other);
    groupElement &operator += (const groupElement &other);
    groupElement &operator *= (const fieldElement &other);

    void print() const;

    bool isOnCurve();
private:
    p224_felem x, y, z;
    p224_felem p_pre_comp[17][3];

    void makePrecomp(p224_felem out[17][3]) const;
};

}
#endif //HYRAX_P224_GROUPELEMENT_HPP
