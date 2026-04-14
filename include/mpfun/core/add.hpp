/**
 * @file add.hpp
 * @brief Addition and subtraction operations for MPFloat
 * 
 * This file implements:
 * - mpnorm: Normalization (carry propagation, rounding, truncation)
 * - mproun: Rounding helper
 * - mpadd: Multi-precision addition
 * - mpsub: Multi-precision subtraction
 * 
 * Ported from mpfunb.f90 subroutines mpadd, mpsub, mpnorm, mproun.
 * 
 * Copyright (c) 2024 Argonne National Laboratory
 * Ported from MPFUN2020-Fort by David H. Bailey
 */

#ifndef MPFUN_CORE_ADD_HPP
#define MPFUN_CORE_ADD_HPP

#include "representation.hpp"

namespace mpfun {
namespace detail {

/**
 * @brief Perform rounding and truncation on an MPFloat
 * 
 * This performs rounding and truncation of the MPR number A. It is called
 * by normalize, and also by other subroutines when the precision level is
 * modified.
 * 
 * Ported from mproun in mpfunb.f90.
 * 
 * @tparam WORDS Precision in mantissa words
 * @param a Array to round (modified in place)
 * @param mpnw Working precision in words
 */
template <int WORDS>
KOKKOS_INLINE_FUNCTION
void round(int64_t* a, int mpnw) {
    int64_t a2 = a[IDX_EXPONENT];
    a[IDX_EXPONENT] = 0;
    
    int ia = extract_sign(a[IDX_SIGN_LENGTH]);
    int na = extract_length(a[IDX_SIGN_LENGTH]);
    if (na > mpnw) na = mpnw;
    int n4 = na + 4;
    
    // Check for initial zeroes
    if (a[IDX_MANTISSA_START] == 0) {
        // Find the first nonzero word and shift left
        int k = 0;
        for (int i = IDX_MANTISSA_START; i <= n4; ++i) {
            if (a[i] != 0) {
                k = i - IDX_MANTISSA_START;
                break;
            }
        }
        
        if (k == 0 && a[IDX_MANTISSA_START] == 0) {
            // All zeros
            a[IDX_PRECISION] = mpnw;
            a[IDX_SIGN_LENGTH] = 0;
            a[IDX_EXPONENT] = 0;
            a[IDX_MANTISSA_START] = 0;
            a[IDX_MANTISSA_START + 1] = 0;
            return;
        }
        
        if (k > 0) {
            // Shift left by k words
            for (int i = IDX_MANTISSA_START; i <= n4 - k; ++i) {
                a[i] = a[i + k];
            }
            a2 = a2 - k;
            na = na - ((k > 2) ? k - 2 : 0);
            if (k == 2) a[na + IDX_MANTISSA_START] = 0;
        }
    }
    
    // Perform rounding
    if (na == mpnw) {
        // Round based on the word beyond precision
        if (a[na + IDX_MANTISSA_START] >= RADIX / 2) {
            a[na + IDX_MANTISSA_START - 1] = a[na + IDX_MANTISSA_START - 1] + 1;
        }
        
        // Release carries due to rounding
        for (int i = na + IDX_MANTISSA_START - 1; i >= IDX_MANTISSA_START; --i) {
            if (a[i] < RADIX) break;
            a[i] = a[i] - RADIX;
            if (i > IDX_MANTISSA_START) {
                a[i - 1] = a[i - 1] + 1;
            } else {
                // Carries propagated all the way - number was all 9's equivalent
                a[IDX_MANTISSA_START] = 1;
                na = 1;
                a2 = a2 + 1;
            }
        }
    }
    
    // Check if trailing word is zero and adjust length
    if (a[na + IDX_MANTISSA_START - 1] == 0) {
        // Find the last nonzero word
        int last_nonzero = IDX_MANTISSA_START - 1;
        for (int i = na + IDX_MANTISSA_START - 1; i >= IDX_MANTISSA_START; --i) {
            if (a[i] != 0) {
                last_nonzero = i;
                break;
            }
        }
        
        if (last_nonzero < IDX_MANTISSA_START) {
            // Result is zero
            a[IDX_PRECISION] = mpnw;
            a[IDX_SIGN_LENGTH] = 0;
            a[IDX_EXPONENT] = 0;
            a[IDX_MANTISSA_START] = 0;
            a[IDX_MANTISSA_START + 1] = 0;
            return;
        }
        
        na = last_nonzero - IDX_MANTISSA_START + 1;
    }
    
    // Check for zero result
    if (a[IDX_MANTISSA_START] == 0) {
        a[IDX_PRECISION] = mpnw;
        a[IDX_SIGN_LENGTH] = 0;
        a[IDX_EXPONENT] = 0;
        a[IDX_MANTISSA_START] = 0;
        a[IDX_MANTISSA_START + 1] = 0;
    } else {
        a[IDX_PRECISION] = mpnw;
        a[IDX_SIGN_LENGTH] = combine_sign_length(ia, na);
        a[IDX_EXPONENT] = a2;
        a[na + IDX_MANTISSA_START] = 0;
        a[na + IDX_MANTISSA_START + 1] = 0;
    }
}

/**
 * @brief Normalize an MPFloat after arithmetic operations
 * 
 * This converts the MP number in array d to the standard normalized form
 * in a. Handles carry propagation, leading zero removal, and rounding.
 * 
 * Ported from mpnorm in mpfunb.f90.
 * 
 * @tparam WORDS Precision in mantissa words
 * @param d Input array (may be modified during normalization)
 * @param a Output array (normalized result)
 * @param mpnw Working precision in words
 */
template <int WORDS>
KOKKOS_INLINE_FUNCTION
void normalize(int64_t* d, int64_t* a, int mpnw) {
    int ia = extract_sign(d[IDX_SIGN_LENGTH]);
    int na = extract_length(d[IDX_SIGN_LENGTH]);
    if (na > mpnw) na = mpnw;
    
    if (na == 0) {
        a[IDX_PRECISION] = mpnw;
        a[IDX_SIGN_LENGTH] = 0;
        a[IDX_EXPONENT] = 0;
        a[IDX_MANTISSA_START] = 0;
        a[IDX_MANTISSA_START + 1] = 0;
        return;
    }
    
    int n4 = na + 4;
    int64_t a2 = d[IDX_EXPONENT];
    d[IDX_EXPONENT] = 0;
    
    // Release carries - may need to repeat if negatives appear
    // Fortran loop: do i = n4, 3, -1; d(i+1) = ...; so it processes d(4) to d(n4+1)
    // In C++ we process d[4] to d[n4], then add carry to d[3]
    bool repeat;
    do {
        repeat = false;
        int64_t t1 = 0;
        
        // Process mantissa words from end to beginning
        // Fortran: do i = n4, 3, -1 with d(i+1) means processing indices 4 to n4+1
        // But we only have mantissa from 4 to na+3 (which is n4-1)
        // Actually the Fortran processes d(4) to d(n4+1) where n4=na+4
        // So the mantissa range is d(4)..d(na+4) in Fortran
        for (int i = n4; i >= IDX_MANTISSA_START; --i) {
            int64_t t3 = t1 + d[i];
            t1 = arith_shift_right(t3, BITS_PER_WORD);
            d[i] = t3 - shift_left(t1, BITS_PER_WORD);
        }
        
        // Add final carry to the exponent overflow slot
        d[IDX_EXPONENT] = d[IDX_EXPONENT] + t1;
        
        if (d[IDX_EXPONENT] < 0) {
            // d[3] is negative -- negate all words and re-normalize
            ia = -ia;
            d[IDX_MANTISSA_START] = d[IDX_MANTISSA_START] + RADIX * d[IDX_EXPONENT];
            d[IDX_EXPONENT] = 0;
            
            // Negate all mantissa words (and sign_length for proper iteration)
            for (int i = IDX_MANTISSA_START; i <= n4; ++i) {
                d[i] = -d[i];
            }
            repeat = true;
        }
    } while (repeat);
    
    if (d[IDX_EXPONENT] > 0) {
        // Nonzero spilled into d[3] - shift right one cell
        for (int i = n4; i >= IDX_EXPONENT; --i) {
            a[i + 1] = d[i];
        }
        a[IDX_MANTISSA_START] = d[IDX_EXPONENT];
        
        na = (na + 1 < mpnw) ? na + 1 : mpnw;
        a2 = a2 + 1;
    } else {
        for (int i = IDX_EXPONENT; i <= n4; ++i) {
            a[i] = d[i];
        }
    }
    
    // Set metadata and round
    a[IDX_ALLOCATED] = mpnw + 6;
    a[IDX_PRECISION] = mpnw;
    a[IDX_SIGN_LENGTH] = combine_sign_length(ia, na);
    a[IDX_EXPONENT] = a2;
    
    round<WORDS>(a, mpnw);
}

/**
 * @brief Add two multi-precision numbers
 * 
 * This routine adds MPR numbers A and B to yield C.
 * 
 * Ported from mpadd in mpfunb.f90.
 * 
 * @tparam WORDS Precision in mantissa words
 * @param a First operand
 * @param b Second operand
 * @param c Result (a + b)
 * @param mpnw Working precision in words
 */
template <int WORDS>
KOKKOS_INLINE_FUNCTION
void add(const int64_t* a, const int64_t* b, int64_t* c, int mpnw) {
    int ia = extract_sign(a[IDX_SIGN_LENGTH]);
    int ib = extract_sign(b[IDX_SIGN_LENGTH]);
    int na = extract_length(a[IDX_SIGN_LENGTH]);
    int nb = extract_length(b[IDX_SIGN_LENGTH]);
    
    if (na > mpnw) na = mpnw;
    if (nb > mpnw) nb = mpnw;
    
    // Check for zero inputs
    if (na == 0) {
        // A is zero -- result is B
        c[IDX_PRECISION] = mpnw;
        c[IDX_SIGN_LENGTH] = combine_sign_length(ib, nb);
        
        for (int i = IDX_SIGN_LENGTH; i <= nb + IDX_EXPONENT; ++i) {
            c[i] = b[i];
        }
        c[nb + IDX_MANTISSA_START] = 0;
        c[nb + IDX_MANTISSA_START + 1] = 0;
        return;
    }
    
    if (nb == 0) {
        // B is zero -- result is A
        c[IDX_PRECISION] = mpnw;
        c[IDX_SIGN_LENGTH] = combine_sign_length(ia, na);
        
        for (int i = IDX_SIGN_LENGTH; i <= na + IDX_EXPONENT; ++i) {
            c[i] = a[i];
        }
        c[na + IDX_MANTISSA_START] = 0;
        c[na + IDX_MANTISSA_START + 1] = 0;
        return;
    }
    
    // Determine if we're adding or subtracting magnitudes
    int idb = (ia == ib) ? 1 : -1;
    
    int64_t ixa = a[IDX_EXPONENT];
    int64_t ixb = b[IDX_EXPONENT];
    int ish = static_cast<int>(ixa - ixb);
    
    // Temporary array for computation
    int64_t d[WORDS + 10];
    for (int i = 0; i < WORDS + 10; ++i) d[i] = 0;
    
    int nd;
    int64_t ixd;
    
    if (ish >= 0) {
        // A has greater or equal exponent -- shift B right
        // m1 = initial A words copied without adding B
        // m2 = end index of A words added to shifted B words
        // m3 = end index of A words copied without adding (after B section)
        // m4 = end index of zero words after A
        // m5 = end index of B words copied with shift
        
        int m1 = (na < ish) ? na : ish;
        int m2 = (na < nb + ish) ? na : nb + ish;
        int m3 = na;
        int m4_temp = (na > ish) ? na : ish;
        int m4 = (m4_temp < mpnw + 1) ? m4_temp : mpnw + 1;
        int m5_temp = (na > nb + ish) ? na : nb + ish;
        int m5 = (m5_temp < mpnw + 1) ? m5_temp : mpnw + 1;
        
        // Copy initial A words
        for (int i = 1; i <= m1; ++i) {
            d[i + IDX_EXPONENT] = a[i + IDX_EXPONENT];
        }
        
        // Add A and shifted B
        for (int i = m1 + 1; i <= m2; ++i) {
            d[i + IDX_EXPONENT] = a[i + IDX_EXPONENT] + idb * b[i - ish + IDX_EXPONENT];
        }
        
        // Copy remaining A words
        for (int i = m2 + 1; i <= m3; ++i) {
            d[i + IDX_EXPONENT] = a[i + IDX_EXPONENT];
        }
        
        // Fill zeros
        for (int i = m3 + 1; i <= m4; ++i) {
            d[i + IDX_EXPONENT] = 0;
        }
        
        // Copy remaining B words
        for (int i = m4 + 1; i <= m5; ++i) {
            d[i + IDX_EXPONENT] = idb * b[i - ish + IDX_EXPONENT];
        }
        
        nd = m5;
        ixd = ixa;
        d[nd + IDX_MANTISSA_START] = 0;
        d[nd + IDX_MANTISSA_START + 1] = 0;
    } else {
        // B has greater exponent -- shift A right
        int nsh = -ish;
        
        int m1 = (nb < nsh) ? nb : nsh;
        int m2 = (nb < na + nsh) ? nb : na + nsh;
        int m3 = nb;
        int m4_temp = (nb > nsh) ? nb : nsh;
        int m4 = (m4_temp < mpnw + 1) ? m4_temp : mpnw + 1;
        int m5_temp = (nb > na + nsh) ? nb : na + nsh;
        int m5 = (m5_temp < mpnw + 1) ? m5_temp : mpnw + 1;
        
        // Copy initial B words (with sign adjustment)
        for (int i = 1; i <= m1; ++i) {
            d[i + IDX_EXPONENT] = idb * b[i + IDX_EXPONENT];
        }
        
        // Add shifted A and B
        for (int i = m1 + 1; i <= m2; ++i) {
            d[i + IDX_EXPONENT] = a[i - nsh + IDX_EXPONENT] + idb * b[i + IDX_EXPONENT];
        }
        
        // Copy remaining B words
        for (int i = m2 + 1; i <= m3; ++i) {
            d[i + IDX_EXPONENT] = idb * b[i + IDX_EXPONENT];
        }
        
        // Fill zeros
        for (int i = m3 + 1; i <= m4; ++i) {
            d[i + IDX_EXPONENT] = 0;
        }
        
        // Copy remaining A words
        for (int i = m4 + 1; i <= m5; ++i) {
            d[i + IDX_EXPONENT] = a[i - nsh + IDX_EXPONENT];
        }
        
        nd = m5;
        ixd = ixb;
        d[nd + IDX_MANTISSA_START] = 0;
        d[nd + IDX_MANTISSA_START + 1] = 0;
    }
    
    // Set up d for normalization and call normalize
    d[IDX_ALLOCATED] = mpnw + 6;
    d[IDX_PRECISION] = mpnw;
    d[IDX_SIGN_LENGTH] = combine_sign_length(ia, nd);
    d[IDX_EXPONENT] = ixd;
    
    normalize<WORDS>(d, c, mpnw);
}

/**
 * @brief Subtract two multi-precision numbers
 * 
 * This routine subtracts MPR numbers A and B to yield C.
 * 
 * Ported from mpsub in mpfunb.f90.
 * 
 * @tparam WORDS Precision in mantissa words
 * @param a First operand
 * @param b Second operand
 * @param c Result (a - b)
 * @param mpnw Working precision in words
 */
template <int WORDS>
KOKKOS_INLINE_FUNCTION
void sub(const int64_t* a, const int64_t* b, int64_t* c, int mpnw) {
    int nb = extract_length(b[IDX_SIGN_LENGTH]);
    if (nb > mpnw) nb = mpnw;
    
    // Create negated copy of b
    int64_t s[WORDS + 10];
    s[IDX_ALLOCATED] = mpnw + 6;
    s[IDX_PRECISION] = mpnw;
    
    if (b[IDX_SIGN_LENGTH] == 0) {
        s[IDX_SIGN_LENGTH] = 0;
    } else if (b[IDX_SIGN_LENGTH] > 0) {
        s[IDX_SIGN_LENGTH] = -nb;
    } else {
        s[IDX_SIGN_LENGTH] = nb;
    }
    
    // Copy exponent and mantissa
    for (int i = IDX_EXPONENT; i <= nb + IDX_MANTISSA_START + 1; ++i) {
        s[i] = b[i];
    }
    
    // Add a + (-b)
    add<WORDS>(a, s, c, mpnw);
}

} // namespace detail
} // namespace mpfun

#endif // MPFUN_CORE_ADD_HPP
