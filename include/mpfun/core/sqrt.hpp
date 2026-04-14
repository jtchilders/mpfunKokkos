/**
 * @file sqrt.hpp
 * @brief Square root operation for MPFloat
 * 
 * This file implements:
 * - mpsqrt: Multi-precision square root using Newton-Raphson iteration
 * 
 * The square root uses the iteration for 1/sqrt(A):
 *    X_{k+1} = X_k + 0.5 * (1 - X_k^2 * A) * X_k
 * 
 * which converges to 1/sqrt(A), then computes sqrt(A) = A * (1/sqrt(A)).
 * 
 * Ported from mpfunb.f90 subroutine mpsqrt.
 * 
 * Copyright (c) 2024 Argonne National Laboratory
 * Ported from MPFUN2020-Fort by David H. Bailey
 */

#ifndef MPFUN_CORE_SQRT_HPP
#define MPFUN_CORE_SQRT_HPP

#include "representation.hpp"
#include "add.hpp"
#include "mul.hpp"
#include "div.hpp"
#include <cmath>

namespace mpfun {
namespace detail {

/**
 * @brief Compute the square root of a multi-precision number
 * 
 * Computes b = sqrt(a) using Newton-Raphson iteration for 1/sqrt(a),
 * then multiplying by a.
 * 
 * The iteration for 1/sqrt(A):
 *    X_{k+1} = X_k + 0.5 * (1 - X_k^2 * A) * X_k
 * 
 * Final iteration uses Karp's trick:
 *    sqrt(A) = (A * X_n) + 0.5 * [A - (A * X_n)^2] * X_n
 * 
 * Ported from mpsqrt in mpfunb.f90.
 * 
 * @tparam WORDS Precision in mantissa words
 * @param a Input (must be non-negative)
 * @param b Output (sqrt(a))
 * @param mpnw Working precision in words
 */
template <int WORDS>
KOKKOS_INLINE_FUNCTION
void sqrt(const int64_t* a, int64_t* b, int mpnw) {
    // Number of iterations before the last one that are repeated
    constexpr int NIT = 3;
    
    int ia = extract_sign(a[IDX_SIGN_LENGTH]);
    int na = extract_length(a[IDX_SIGN_LENGTH]);
    
    if (na > mpnw) na = mpnw;
    
    // Check for zero input
    if (na == 0) {
        b[IDX_PRECISION] = mpnw;
        b[IDX_SIGN_LENGTH] = 0;
        b[IDX_EXPONENT] = 0;
        b[IDX_MANTISSA_START] = 0;
        b[IDX_MANTISSA_START + 1] = 0;
        return;
    }
    
    // Check for negative input - error
    if (ia < 0) {
        b[IDX_ALLOCATED] = -1;  // NaN marker
        b[IDX_PRECISION] = mpnw;
        b[IDX_SIGN_LENGTH] = 0;
        return;
    }
    
    // Temporary arrays
    int64_t s0[WORDS + 10];
    int64_t s1[WORDS + 10];
    int64_t s2[WORDS + 10];
    int64_t s3[WORDS + 10];
    
    for (int i = 0; i < WORDS + 10; ++i) {
        s0[i] = s1[i] = s2[i] = s3[i] = 0;
    }
    
    s0[IDX_ALLOCATED] = mpnw + 7;
    s1[IDX_ALLOCATED] = mpnw + 7;
    s2[IDX_ALLOCATED] = mpnw + 7;
    s3[IDX_ALLOCATED] = mpnw + 7;
    
    // Determine the least integer MQ such that 2^MQ >= MPNW
    double t1 = static_cast<double>(mpnw);
    int mq = static_cast<int>(LOG2_E * Kokkos::log(t1) + 1.0 - COMPARE_FUZZ);
    
    // Compute the initial approximation of 1 / sqrt(A)
    double da;
    int n_exp;
    mpmdc<WORDS>(a, da, n_exp, mpnw);
    
    // Compute 1/sqrt(da * 2^n)
    // n2 = -n/2, then adjust da to account for odd n
    int n2 = -n_exp / 2;
    double t2 = Kokkos::sqrt(da * Kokkos::pow(2.0, static_cast<double>(n_exp + 2 * n2)));
    t1 = 1.0 / t2;
    mpdmc<WORDS>(t1, n2, s2, mpnw);
    
    // Initialize s3 = 1.0
    mpdmc<WORDS>(1.0, 0, s3, mpnw);
    
    int mpnw1 = 5;
    int iq = 0;
    int nw1 = mpnw1;
    int nw2 = mpnw1;
    
    // Newton-Raphson iteration: X_{k+1} = X_k + 0.5 * (1 - X_k^2 * A) * X_k
    for (int k = 1; k <= mq - 1; ++k) {
        if (k > 2) {
            nw1 = mpnw1;
            mpnw1 = (2 * mpnw1 - 2 < mpnw) ? 2 * mpnw1 - 2 : mpnw;
            mpnw1 += 1;
            nw2 = mpnw1;
        }
        
        do {
            mul<WORDS>(s2, s2, s0, nw2);     // s0 = X_k^2
            mul<WORDS>(a, s0, s1, nw2);      // s1 = A * X_k^2
            sub<WORDS>(s3, s1, s0, nw2);     // s0 = 1 - A * X_k^2
            mul<WORDS>(s2, s0, s1, nw1);     // s1 = X_k * (1 - A * X_k^2)
            muld<WORDS>(s1, 0.5, s0, nw1);   // s0 = 0.5 * X_k * (1 - A * X_k^2)
            add<WORDS>(s2, s0, s1, nw2);     // s1 = X_k + 0.5 * X_k * (1 - A * X_k^2)
            mpeq<WORDS>(s1, s2, nw2);        // s2 = s1
            
            if (k == mq - NIT && iq == 0) {
                iq = 1;
            } else {
                break;
            }
        } while (true);
    }
    
    // Final iteration using Karp's trick:
    // sqrt(A) = (A * X_n) + 0.5 * [A - (A * X_n)^2] * X_n
    
    nw1 = mpnw1;
    mpnw1 = (2 * mpnw1 - 2 < mpnw) ? 2 * mpnw1 - 2 : mpnw;
    mpnw1 += 1;
    nw2 = mpnw1;
    
    mul<WORDS>(a, s2, s0, nw1);       // s0 = A * X_n
    mul<WORDS>(s0, s0, s1, nw2);      // s1 = (A * X_n)^2
    sub<WORDS>(a, s1, s3, nw2);       // s3 = A - (A * X_n)^2
    mul<WORDS>(s3, s2, s1, nw1);      // s1 = [A - (A * X_n)^2] * X_n
    muld<WORDS>(s1, 0.5, s3, nw1);    // s3 = 0.5 * [A - (A * X_n)^2] * X_n
    add<WORDS>(s0, s3, s2, nw2);      // s2 = (A * X_n) + 0.5 * [A - (A * X_n)^2] * X_n
    
    // Round and copy to result
    round<WORDS>(s2, mpnw);
    mpeq<WORDS>(s2, b, mpnw);
}

} // namespace detail

// Free function for MPFloat class
template <int WORDS>
class MPFloat;

/**
 * @brief Compute the square root of an MPFloat
 * 
 * @tparam W Precision in mantissa words
 * @param x Input value (must be non-negative)
 * @return Square root of x
 */
template <int W>
KOKKOS_INLINE_FUNCTION
MPFloat<W> sqrt(const MPFloat<W>& x) {
    MPFloat<W> result;
    detail::sqrt<W>(x.data(), result.data(), W);
    return result;
}

} // namespace mpfun

#endif // MPFUN_CORE_SQRT_HPP
