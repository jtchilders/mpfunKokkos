/**
 * @file div.hpp
 * @brief Division operations for MPFloat
 * 
 * This file implements:
 * - mpdiv: Multi-precision division using Newton-Raphson iteration
 * - mpdivd: Multi-precision divided by double
 * 
 * The division uses the iteration:
 *    X_{k+1} = X_k + (1 - X_k * B) * X_k
 * 
 * which converges to 1/B, followed by multiplication by A.
 * 
 * Ported from mpfunb.f90 subroutine mpdiv.
 * 
 * Copyright (c) 2024 Argonne National Laboratory
 * Ported from MPFUN2020-Fort by David H. Bailey
 */

#ifndef MPFUN_CORE_DIV_HPP
#define MPFUN_CORE_DIV_HPP

#include "representation.hpp"
#include "add.hpp"
#include "mul.hpp"
#include <cmath>

namespace mpfun {
namespace detail {

/**
 * @brief Extract a double approximation from an MPFloat
 * 
 * Returns b and n such that the MPFloat value ≈ b * 2^n
 * with 1 <= |b| < 2.
 * 
 * Ported from mpmdc in mpfunb.f90.
 * 
 * @tparam WORDS Precision in mantissa words
 * @param a Input MPFloat array
 * @param b Output double approximation (in range [1,2) or (-2,-1])
 * @param n Output power of 2 exponent
 * @param mpnw Working precision
 */
template <int WORDS>
KOKKOS_INLINE_FUNCTION
void mpmdc(const int64_t* a, double& b, int& n, int mpnw) {
    if (a[IDX_SIGN_LENGTH] == 0) {
        b = 0.0;
        n = 0;
        return;
    }
    
    int na = extract_length(a[IDX_SIGN_LENGTH]);
    int sa = extract_sign(a[IDX_SIGN_LENGTH]);
    
    // Get double approximation from mantissa
    double aa = static_cast<double>(a[IDX_MANTISSA_START]);
    if (na >= 2) {
        aa += static_cast<double>(a[IDX_MANTISSA_START + 1]) / static_cast<double>(RADIX);
    }
    
    // Compute power of 2 from exponent
    n = static_cast<int>(BITS_PER_WORD * a[IDX_EXPONENT]);
    b = (sa >= 0) ? aa : -aa;
    
    // Reduce b to within [1, 2)
    int na2 = static_cast<int>(Kokkos::log(Kokkos::fabs(b)) / Kokkos::log(2.0) + COMPARE_FUZZ);
    b = b / Kokkos::pow(2.0, static_cast<double>(na2));
    n = n + na2;
    
    if (Kokkos::fabs(b) < 1.0) {
        b = 2.0 * b;
        n = n - 1;
    } else if (Kokkos::fabs(b) >= 2.0) {
        b = 0.5 * b;
        n = n + 1;
    }
}

/**
 * @brief Convert double * 2^n to MPFloat format
 * 
 * Converts a * 2^n to MPFloat representation.
 * 
 * Ported from mpdmc in mpfunb.f90.
 * 
 * @tparam WORDS Precision in mantissa words
 * @param a Double precision value
 * @param n Power of 2 scaling
 * @param b Output MPFloat array
 * @param mpnw Working precision
 */
template <int WORDS>
KOKKOS_INLINE_FUNCTION
void mpdmc(double a, int n, int64_t* b, int mpnw) {
    b[IDX_ALLOCATED] = mpnw + 6;
    b[IDX_PRECISION] = mpnw;
    
    if (a == 0.0) {
        b[IDX_SIGN_LENGTH] = 0;
        b[IDX_EXPONENT] = 0;
        b[IDX_MANTISSA_START] = 0;
        b[IDX_MANTISSA_START + 1] = 0;
        return;
    }
    
    int n1 = n / BITS_PER_WORD;
    int n2 = n - BITS_PER_WORD * n1;
    double aa = Kokkos::fabs(a) * Kokkos::pow(2.0, static_cast<double>(n2));
    
    const double bdx = static_cast<double>(RADIX);
    
    // Reduce aa to within [1, RADIX)
    if (aa >= bdx) {
        for (int k = 1; k <= 100; ++k) {
            aa = aa / bdx;
            if (aa < bdx) {
                n1 = n1 + k;
                break;
            }
        }
    } else if (aa < 1.0) {
        for (int k = 1; k <= 100; ++k) {
            aa = aa * bdx;
            if (aa >= 1.0) {
                n1 = n1 - k;
                break;
            }
        }
    }
    
    // Store result
    b[IDX_EXPONENT] = n1;
    b[IDX_MANTISSA_START] = static_cast<int64_t>(aa);
    aa = bdx * (aa - static_cast<double>(b[IDX_MANTISSA_START]));
    b[IDX_MANTISSA_START + 1] = static_cast<int64_t>(aa);
    b[IDX_MANTISSA_START + 2] = 0;
    b[IDX_MANTISSA_START + 3] = 0;
    b[IDX_MANTISSA_START + 4] = 0;
    
    // Determine actual length (2 words from double)
    int len = 1;
    if (b[IDX_MANTISSA_START + 1] != 0) len = 2;
    
    b[IDX_SIGN_LENGTH] = (a >= 0.0) ? len : -len;
}

/**
 * @brief Copy one MPFloat to another
 * 
 * Ported from mpeq in mpfunb.f90.
 * 
 * @tparam WORDS Precision in mantissa words
 * @param a Source array
 * @param b Destination array
 * @param mpnw Working precision
 */
template <int WORDS>
KOKKOS_INLINE_FUNCTION
void mpeq(const int64_t* a, int64_t* b, int mpnw) {
    int na = extract_length(a[IDX_SIGN_LENGTH]);
    if (na > mpnw) na = mpnw;
    
    b[IDX_PRECISION] = mpnw;
    b[IDX_SIGN_LENGTH] = a[IDX_SIGN_LENGTH];
    
    for (int i = IDX_EXPONENT; i <= na + IDX_EXPONENT; ++i) {
        b[i] = a[i];
    }
    
    b[na + IDX_MANTISSA_START] = 0;
    b[na + IDX_MANTISSA_START + 1] = 0;
}

/**
 * @brief Divide two multi-precision numbers
 * 
 * Computes c = a / b using Newton-Raphson iteration for 1/b,
 * then multiplying by a.
 * 
 * Ported from mpdiv in mpfunb.f90.
 * 
 * @tparam WORDS Precision in mantissa words
 * @param a Dividend
 * @param b Divisor
 * @param c Result (a / b)
 * @param mpnw Working precision in words
 */
template <int WORDS>
KOKKOS_INLINE_FUNCTION
void div(const int64_t* a, const int64_t* b, int64_t* c, int mpnw) {
    // Number of iterations before the last one that are repeated
    constexpr int NIT = 3;
    
    int na = extract_length(a[IDX_SIGN_LENGTH]);
    int nb = extract_length(b[IDX_SIGN_LENGTH]);
    
    if (na > mpnw) na = mpnw;
    if (nb > mpnw) nb = mpnw;
    
    // Check for zero dividend
    if (na == 0) {
        c[IDX_PRECISION] = mpnw;
        c[IDX_SIGN_LENGTH] = 0;
        c[IDX_EXPONENT] = 0;
        c[IDX_MANTISSA_START] = 0;
        c[IDX_MANTISSA_START + 1] = 0;
        return;
    }
    
    // Check for zero divisor - this is an error
    if (nb == 0) {
        c[IDX_ALLOCATED] = -1;  // NaN marker
        c[IDX_PRECISION] = mpnw;
        c[IDX_SIGN_LENGTH] = 0;
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
    
    // Compute the initial approximation of 1 / B using mpmdc/mpdmc pattern
    double db;
    int n_exp;
    mpmdc<WORDS>(b, db, n_exp, mpnw);
    
    // Compute 1/db and convert back to MP format
    double t2 = 1.0 / db;
    mpdmc<WORDS>(t2, -n_exp, s2, mpnw);
    
    // Initialize s3 = 1.0
    mpdmc<WORDS>(1.0, 0, s3, mpnw);
    
    int mpnw1 = 5;
    int iq = 0;
    int nw1 = mpnw1;
    int nw2 = mpnw1;
    
    // Newton-Raphson iteration: X_{k+1} = X_k + (1 - X_k * B) * X_k
    for (int k = 1; k <= mq - 1; ++k) {
        if (k > 2) {
            nw1 = mpnw1;
            mpnw1 = (2 * mpnw1 - 2 < mpnw) ? 2 * mpnw1 - 2 : mpnw;
            mpnw1 += 1;
            nw2 = mpnw1;
        }
        
        do {
            mul<WORDS>(b, s2, s1, nw2);      // s1 = B * X_k
            sub<WORDS>(s3, s1, s0, nw2);     // s0 = 1 - B * X_k
            mul<WORDS>(s2, s0, s1, nw1);     // s1 = X_k * (1 - B * X_k)
            add<WORDS>(s2, s1, s0, nw2);     // s0 = X_k + X_k * (1 - B * X_k)
            mpeq<WORDS>(s0, s2, nw2);        // s2 = s0
            
            if (k == mq - NIT && iq == 0) {
                iq = 1;
            } else {
                break;
            }
        } while (true);
    }
    
    // Final iteration using Karp's trick:
    // A/B = (A * X_n) + [A - (A * X_n) * B] * X_n
    
    nw1 = mpnw1;
    mpnw1 = (2 * mpnw1 - 1 < mpnw) ? 2 * mpnw1 - 1 : mpnw;
    mpnw1 += 1;
    nw2 = mpnw1;
    
    mul<WORDS>(a, s2, s0, nw1);      // s0 = A * X_n
    mul<WORDS>(s0, b, s1, nw2);      // s1 = (A * X_n) * B
    sub<WORDS>(a, s1, s3, nw2);      // s3 = A - (A * X_n) * B
    mul<WORDS>(s3, s2, s1, nw1);     // s1 = [A - (A * X_n) * B] * X_n
    add<WORDS>(s0, s1, s2, nw2);     // s2 = (A * X_n) + [A - (A * X_n) * B] * X_n
    
    // Round and copy to result
    round<WORDS>(s2, mpnw);
    mpeq<WORDS>(s2, c, mpnw);
}

/**
 * @brief Divide multi-precision number by double
 * 
 * Computes c = a / b where a is MPFloat and b is double.
 * 
 * @tparam WORDS Precision in mantissa words
 * @param a MPFloat dividend
 * @param b Double divisor
 * @param c Result
 * @param mpnw Working precision
 */
template <int WORDS>
KOKKOS_INLINE_FUNCTION
void divd(const int64_t* a, double b, int64_t* c, int mpnw) {
    // Convert b to MPFloat and use full division
    int64_t bd[WORDS + 10];
    bd[IDX_ALLOCATED] = mpnw + 6;
    
    if (b == 0.0) {
        // Division by zero - set NaN
        c[IDX_ALLOCATED] = -1;
        c[IDX_PRECISION] = mpnw;
        c[IDX_SIGN_LENGTH] = 0;
        return;
    }
    
    mpdmc<WORDS>(b, 0, bd, mpnw);
    div<WORDS>(a, bd, c, mpnw);
}

} // namespace detail
} // namespace mpfun

#endif // MPFUN_CORE_DIV_HPP
