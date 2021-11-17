# hyrax-bls12-381

## Introduction
This is a polynomial commitment implementation refer to [Hyrax](https://eprint.iacr.org/2017/1132.pdf) based on BLS12-381 implemented by [mcl](https://github.com/herumi/mcl). The underline operations of scalar and group element refers to OpenSSL.
This scheme is particularly for multi-linear extension of an array.

## Usage
You can use this polynomial commitment in this way:
```C++

#include <bits/stdc++.h>

// Include the header
#include <hyrax-bls12-381/polyCommit.hpp>

#define F Fr  // This is the field on which the polynomial is defined.
#define G G1  // BLS12-381 is friendly to pairing, but here we only use G1.
#define F_ONE (Fr::one())
#define F_ZERO (Fr(0))

int main() {
    u64 n = 1 << 10;
    int logn = 10;
    
    // Initialize the polynomial coefficients and evaluation point.
    vector<F> polynomial(n), eval_point(logn);
    for (auto &x: polynomial) x.setByCSPRNG();
    for (auto &x: eval_point) x.setByCSPRNG();
    
    // Initialize the structure setting.
    initPairing(mcl::BLS12_381);

    // Generate the generators used in Hyrax.
    u64 n_sqrt = 1ULL << (logn - (logn >> 1));
    vector<G> gens(n_sqrt);
    for (auto &x: gens) {
        F tmp;
        tmp.setByCSPRNG();
        x = mcl::bn::getG1basePoint() * tmp;
    }
    
    // Initialize the prover.
    hyrax_bls12_381::polyProver poly_p(polynomial, gens);
    
    // Calculate the evaluation
    F eval_value = poly_p.evaluate(eval_point);
    
    // Initialize the verifier and verify the evaluation value.
    hyrax_bls12_381::polyVerifier poly_v(poly_p, gens);
    if (poly_v.verify(eval_point, eval_value)) std::cout << "Correct!" << std::endl;
    else std::cout << "Incorrect!" << std::endl;
    
    return 0;
}
```

## Reference
- [Doubly-efficient zksnarks without trusted setup](https://doi.org/10.1109/SP.2018.00060). Wahby, R. S., Tzialla, I., Shelat, A., Thaler, J., & Walfish, M. (S&P 2018)

- [Hyrax](https://github.com/hyraxZK/hyraxZK.git)

- [mcl](https://github.com/herumi/mcl)
