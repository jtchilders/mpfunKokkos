/**
 * @file mul.hpp
 * @brief Multiplication operations for MPFloat
 * 
 * This file implements:
 * - mpmul: Multi-precision multiplication (schoolbook algorithm)
 * - mpmuld: Multi-precision times double
 * 
 * TODO: Add FFT-based multiplication for very high precision (>20k digits)
 * 
 * Ported from mpfunb.f90 subroutine mpmul.
 * 
 * Copyright (c) 2024 Argonne National Laboratory
 * Ported from MPFUN2020-Fort by David H. Bailey
 */

#ifndef MPFUN_CORE_MUL_HPP
#define MPFUN_CORE_MUL_HPP

#include "representation.hpp"
#include "add.hpp"

namespace mpfun {
namespace detail {

/**
 * @brief Multiply two multi-precision numbers
 * 
 * This routine multiplies MPR numbers A and B to yield C.
 * Uses the schoolbook O(n²) algorithm, which is efficient for
 * moderate precision (up to a few thousand digits).
 * 
 * Ported from mpmul in mpfunb.f90.
 * 
 * @tparam WORDS Precision in mantissa words
 * @param a First operand
 * @param b Second operand
 * @param c Result (a * b)
 * @param mpnw Working precision in words
 */
template <int WORDS>
KOKKOS_INLINE_FUNCTION
void mul(const int64_t* a, const int64_t* b, int64_t* c, int mpnw) {
    int ia = extract_sign(a[IDX_SIGN_LENGTH]);
    int ib = extract_sign(b[IDX_SIGN_LENGTH]);
    int na = extract_length(a[IDX_SIGN_LENGTH]);
    int nb = extract_length(b[IDX_SIGN_LENGTH]);
    
    if (na > mpnw) na = mpnw;
    if (nb > mpnw) nb = mpnw;
    
    int nc = na + nb;
    if (nc > mpnw) nc = mpnw;
    
    // Temporary array for result
    int64_t d[2 * WORDS + 10];
    for (int i = 0; i < 2 * WORDS + 10; ++i) d[i] = 0;
    
    // Check for zero inputs
    if (na == 0 || nb == 0) {
        c[IDX_PRECISION] = mpnw;
        c[IDX_SIGN_LENGTH] = 0;
        c[IDX_EXPONENT] = 0;
        c[IDX_MANTISSA_START] = 0;
        c[IDX_MANTISSA_START + 1] = 0;
        return;
    }
    
    // Check for multiply by 1 or -1
    if (na == 1 && a[IDX_MANTISSA_START] == 1) {
        // A is ±1
        c[IDX_PRECISION] = mpnw;
        c[IDX_SIGN_LENGTH] = combine_sign_length(ia * ib, nb);
        c[IDX_EXPONENT] = a[IDX_EXPONENT] + b[IDX_EXPONENT];
        for (int i = 3; i <= nb + 2; ++i) {
            c[i + 1] = b[i + 1];
        }
        c[nb + 4] = 0;
        c[nb + 5] = 0;
        return;
    }
    
    if (nb == 1 && b[IDX_MANTISSA_START] == 1) {
        // B is ±1
        c[IDX_PRECISION] = mpnw;
        c[IDX_SIGN_LENGTH] = combine_sign_length(ia * ib, na);
        c[IDX_EXPONENT] = a[IDX_EXPONENT] + b[IDX_EXPONENT];
        for (int i = 3; i <= na + 2; ++i) {
            c[i + 1] = a[i + 1];
        }
        c[na + 4] = 0;
        c[na + 5] = 0;
        return;
    }
    
    // Initialize result
    int64_t dd = a[IDX_EXPONENT] + b[IDX_EXPONENT];
    d[IDX_ALLOCATED] = mpnw + 6;
    d[IDX_PRECISION] = mpnw;
    d[IDX_SIGN_LENGTH] = combine_sign_length(ia * ib, nc);
    
    // Clear result mantissa
    for (int i = 2; i <= nc + 5; ++i) {
        d[i + 1] = 0;
    }
    
    // Schoolbook multiplication with half-word splitting
    // This avoids 64-bit overflow during accumulation
    
    for (int j = 3; j <= na + 2; ++j) {
        int j3 = j - 3;
        int n2 = (nb + 2 < mpnw + 4 - j3) ? nb + 2 : mpnw + 4 - j3;
        
        int64_t aj = a[j + 1];
        int64_t a1, a2;
        split_word(aj, a1, a2);
        
        for (int i = 3; i <= n2; ++i) {
            int64_t bi = b[i + 1];
            int64_t b1, b2;
            split_word(bi, b1, b2);
            
            // Compute product using Karatsuba-like decomposition
            // aj * bi = (a1*2^30 + a2) * (b1*2^30 + b2)
            //         = a1*b1*2^60 + (a1*b2 + a2*b1)*2^30 + a2*b2
            
            int64_t c1 = a1 * b2 + a2 * b1;
            int64_t c2 = arith_shift_right(c1, HALF_BITS);
            int64_t c3 = c1 - shift_left(c2, HALF_BITS);
            
            d[i + j3] = d[i + j3] + a1 * b1 + c2;
            d[i + j3 + 1] = d[i + j3 + 1] + a2 * b2 + shift_left(c3, HALF_BITS);
        }
        
        // Release carries on the just-computed section
        int64_t t1 = 0;
        for (int i = n2; i >= 3; --i) {
            int64_t t3 = t1 + d[i + j3 + 1];
            t1 = arith_shift_right(t3, BITS_PER_WORD);
            d[i + j3 + 1] = t3 - shift_left(t1, BITS_PER_WORD);
        }
        d[j3 + 3] = d[j3 + 3] + t1;
    }
    
    // Release carries on the full d vector
    int64_t t1 = 0;
    for (int i = nc + 1; i >= 1; --i) {
        int64_t t3 = t1 + d[i + 3];
        t1 = arith_shift_right(t3, BITS_PER_WORD);
        d[i + 3] = t3 - shift_left(t1, BITS_PER_WORD);
    }
    d[3] = d[3] + t1;
    
    // If d[3] is nonzero, shift right one cell
    if (d[3] != 0) {
        dd = dd + 1;
        nc = (nc + 1 < mpnw) ? nc + 1 : mpnw;
        d[IDX_SIGN_LENGTH] = combine_sign_length(ia * ib, nc);
        
        for (int i = nc + 4; i >= 3; --i) {
            d[i + 1] = d[i];
        }
    }
    
    d[IDX_EXPONENT] = dd;
    
    // Copy and normalize
    for (int i = 1; i <= nc + 5; ++i) {
        c[i] = d[i];
    }
    
    normalize<WORDS>(c, c, mpnw);
}

/**
 * @brief Multiply multi-precision number by double
 * 
 * Computes c = a * b where a is MPFloat and b is double.
 * 
 * Note: The product is fully accurate only if b is an exact binary value.
 * 
 * @tparam WORDS Precision in mantissa words
 * @param a MPFloat operand
 * @param b Double operand
 * @param c Result
 * @param mpnw Working precision
 */
template <int WORDS>
KOKKOS_INLINE_FUNCTION
void muld(const int64_t* a, double b, int64_t* c, int mpnw) {
    // For now, convert b to MPFloat and use full multiplication
    // TODO: Optimize this using short multiplication
    int64_t bd[WORDS + 10];
    bd[IDX_ALLOCATED] = mpnw + 6;
    bd[IDX_PRECISION] = mpnw;
    
    // Simple conversion of b to MPFloat format
    if (b == 0.0) {
        c[IDX_PRECISION] = mpnw;
        c[IDX_SIGN_LENGTH] = 0;
        c[IDX_EXPONENT] = 0;
        c[IDX_MANTISSA_START] = 0;
        return;
    }
    
    double bb = Kokkos::fabs(b);
    int n1 = 0;
    const double bdx = static_cast<double>(RADIX);
    
    // Reduce bb to [1, RADIX)
    while (bb >= bdx) {
        bb /= bdx;
        ++n1;
    }
    while (bb < 1.0) {
        bb *= bdx;
        --n1;
    }
    
    bd[IDX_SIGN_LENGTH] = (b >= 0.0) ? 2 : -2;
    bd[IDX_EXPONENT] = n1;
    bd[IDX_MANTISSA_START] = static_cast<int64_t>(bb);
    bb = bdx * (bb - static_cast<double>(bd[IDX_MANTISSA_START]));
    bd[IDX_MANTISSA_START + 1] = static_cast<int64_t>(bb);
    bd[IDX_MANTISSA_START + 2] = 0;
    bd[IDX_MANTISSA_START + 3] = 0;
    
    mul<WORDS>(a, bd, c, mpnw);
}

} // namespace detail
} // namespace mpfun

#endif // MPFUN_CORE_MUL_HPP
