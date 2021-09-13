//
// Created by juzix on 2021/6/1.
//

#include "fieldElement.hpp"

using std::move;
using std::is_same;

namespace hyrax_p224 {
    const u8 mx_ele[28] = {0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 127, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 127};
    const u8 mod[28]    = {0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x01};
    const u8 root_of_unity[28] = {220, 212, 90, 196, 98, 37, 229, 179, 141, 59, 55, 47, 9, 62, 249, 238, 114, 177, 35, 123, 237, 233, 18, 112, 66, 69, 140, 241};
    const int ord = 64;

    bool fieldElement::initialized = false;
    int fieldElement::multCounter, fieldElement::addCounter;
    bool fieldElement::isCounting;
    bool fieldElement::isSumchecking;

    void fieldElement::init() {
        initialized = true;
        srand(time(NULL));
        isCounting = false;
    }

    fieldElement fieldElement::get_root_of_unity(int order) {
        fieldElement res(root_of_unity);

        for (int i = ord; i > order; --i)
            res.setSqr();
        return res;
    }

    fieldElement fieldElement::zero() {
        return fieldElement();
    }

    fieldElement fieldElement::one() {
        p224_felem dat;
        p224_felem_one(dat);
        return fieldElement(dat);
    }

    fieldElement fieldElement::random() {
        u8 tmp[28];
        for (int i = 0; i < 28; ++i) tmp[i] = myrandom();
        fieldElement res(tmp);
        p224_widefelem tmp2{res.data[0], res.data[1], res.data[2], res.data[3], 0, 0, 0};
        p224_felem_reduce(res.data, tmp2);
        p224_felem_contract(res.data, res.data);
        return res;
    }

    fieldElement fieldElement::operator+(const fieldElement &other) const {
        if (isCounting) ++addCounter;
        fieldElement res;
        p224_felem_assign(res.data, data);
        p224_felem_mysum(res.data, other.data);
        p224_felem_contract(res.data, res.data);
        return res;
    }

    fieldElement fieldElement::operator-(const fieldElement &other) const {
        if (isCounting) ++addCounter;
        fieldElement res;
        p224_widefelem tmp = {data[0], data[1], data[2], data[3], 0, 0, 0};
        p224_felem_diff_128_64(tmp, other.data);
        p224_felem_reduce(res.data, tmp);
        p224_felem_contract(res.data, res.data);
        return res;
    }

    fieldElement fieldElement::operator-() const {
        if (isCounting) ++addCounter;
        fieldElement res;
        p224_felem_neg(res.data, data);
        p224_felem_contract(res.data, res.data);
        return res;
    }

    bool fieldElement::isNegative() const {
        for (int i = 3; i >= 0; --i) if (mx_ele[i] != data[i]) return mx_ele[i] < data[i];
        return false;
    }

    fieldElement fieldElement::operator*(const fieldElement &other) const {
        if (isCounting) ++multCounter;
        fieldElement res;
        p224_widefelem tmp;

        p224_felem_mul(tmp, data, other.data);
        p224_felem_reduce(res.data, tmp);
        p224_felem_contract(res.data, res.data);
        return res;
    }

    bool fieldElement::operator==(const fieldElement &other) const {
        return p224_felem_equal(data, other.data);
    }

    bool fieldElement::operator!=(const fieldElement &other) const {
        return !p224_felem_equal(data, other.data);
    }

    fieldElement &fieldElement::operator=(const fieldElement &other) {
        p224_felem_assign(data, other.data);
        return *this;
    }

    fieldElement &fieldElement::operator+=(const fieldElement &other) {
        *this = *this + other;
        return *this;
    }

    fieldElement &fieldElement::operator-=(const fieldElement &other) {
        *this = *this - other;
        return *this;
    }

    fieldElement &fieldElement::operator*=(const fieldElement &other) {
        *this = *this * other;
        return *this;
    }

    ostream &operator << (ostream &out, const fieldElement &data) {
        fieldElement tmp = data.isNegative() ? -data : data;
        out << (data.isNegative() ? "-" : "") << "(" << tmp[0] << ' ' << tmp[1]<< ' ' << tmp[2] << ' ' << tmp[3] << ")";
        return out;
    }

    fieldElement fieldElement::innerProd(vector<fieldElement>::iterator a, vector<fieldElement>::iterator b, u64 n) {
        fieldElement res;
        for (int i = 0; i < n; ++i)
            res += (a[i] * b[i]);
        return res;
    }

    u8 fieldElement::getBit(unsigned int i) const {
        return p224_get_bit_from_felem(data, i);
    }

    fieldElement::fieldElement(const u8 *in) {
        p224_bin28_to_felem(data, in);
    }

    fieldElement::fieldElement(const p224_felem in) {
        memcpy(data, in, sizeof data);
    }

    fieldElement fieldElement::maxWithZero(const fieldElement &a, const fieldElement &b) {
        auto res = a.lessThan(b) ? b : a;
        if (res.isNegative()) res = zero();
        return res;
    }
    fieldElement fieldElement::maxUnsigned(const fieldElement &a, const fieldElement &b) {
        auto res = a.lessThan(b) ? b : a;
        return res;
    }

    fieldElement *fieldElement::generateRandomness(unsigned int size) {
        fieldElement *tmp = new fieldElement[size];
        for (int i = 0; i < size; ++i) tmp[i] = random();
        return tmp;
    }

    fieldElement fieldElement::abs() const {
        return isNegative() ? -*this : *this;
    }

    fieldElement fieldElement::sqr() const {
        fieldElement res;
        p224_widefelem tmp;
        p224_felem_square(tmp, data);
        p224_felem_reduce(res.data, tmp);
        p224_felem_contract(res.data, res.data);
        return res;
    }

    fieldElement fieldElement::inv() const {
        fieldElement res;
        p224_felem_inv(res.data, data);
        p224_felem_contract(res.data, res.data);
        return res;
    }

    void fieldElement::setAbs() {
        *this = this -> abs();
    }

    void fieldElement::setSqr() {
        p224_widefelem tmp;
        p224_felem_square(tmp, data);
        p224_felem_reduce(data, tmp);
        p224_felem_contract(data, data);
    }

    void fieldElement::setInv() {
        p224_felem_inv(data, data);
        p224_felem_contract(data, data);
    }

    p224_limb fieldElement::operator[](const u8 i) const {
        return data[i];
    }

    p224_limb &fieldElement::operator[](const u8 i) {
        return data[i];
    }

    fieldElement::operator bool() const {
        return data[0] || data[1] || data[2] || data[3];
    }

    __int128_t fieldElement::toint128() const {
        bool flag = isNegative();
        fieldElement tmp = flag ? -*this : *this;
        auto res = ((u128) tmp.data[1]) << 56 | tmp.data[0];
        return flag ? -res : res;
    }

    bool fieldElement::lessThan(const fieldElement &b) const {
        return (*this - b).isNegative();
    }

    void fieldElement::toBin28(u8 *bytes) const {
        p224_felem_to_bin28(bytes, data);
    }

    void fieldElement::print(FILE *fileno) const {
        fprintf(fileno, "0x%llx", data[3]);
        for (int i = 2; i >= 0; --i) fprintf(fileno, "%014llx", data[i]);
        fprintf(fileno, "\n");
    }

    u8 fieldElement::getBitWidth() const {
        u8 res = data[3] ? 56 * 3 : (data[2] ? 56 * 2 : (data[1] ? 56 : 0));
        u64 dat = data[3] ? data[3] : (data[2] ? data[2] : (data[1] ? data[1] : data[0]));
        if (!dat) return res;
        for (int i = 32; i && dat; i >>= 1) {
            if (dat >> i) {
                res += i;
                dat >>= i;
            }
        }
        return res + 1;
    }

    bool fieldElement::isZero() {
        return p224_felem_is_zero(data);
    }

    char *fieldElement::toString() const {
        std::string res = "0x";
        for (int i = 3; i >= 0; --i) {
            char s[30];
            sprintf(s, "%014llx", data[i]);
            res += std::string(s);
        }
        return strdup(res.c_str());
    }
}