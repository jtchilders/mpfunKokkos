/**
 * @file compare.hpp
 * @brief Comparison operations for MPFloat
 * 
 * This file implements:
 * - mpcpr: Multi-precision comparison
 * 
 * Ported from mpfunb.f90 subroutine mpcpr.
 * 
 * Copyright (c) 2024 Argonne National Laboratory
 * Ported from MPFUN2020-Fort by David H. Bailey
 */

#ifndef MPFUN_CORE_COMPARE_HPP
#define MPFUN_CORE_COMPARE_HPP

#include "representation.hpp"
#include "add.hpp"

namespace mpfun {

// Forward declaration
template <int WORDS> class MPFloat;

/**
 * @brief Compare two multi-precision numbers
 * 
 * This routine compares the MPR numbers A and B and returns:
 *   -1 if A < B
 *    0 if A = B
 *   +1 if A > B
 * 
 * Note that the metadata words do NOT need to be the same for the
 * result to be "equal" - only the mathematical values are compared.
 * 
 * Ported from mpcpr in mpfunb.f90.
 * 
 * @tparam W Precision in mantissa words
 * @param a First operand
 * @param b Second operand
 * @return -1, 0, or +1 indicating comparison result
 */
template <int W>
KOKKOS_INLINE_FUNCTION
int compare(const MPFloat<W>& a, const MPFloat<W>& b) {
    // Quick sign-based comparisons
    int sa = a.sign();
    int sb = b.sign();
    
    // If signs differ, comparison is immediate
    if (sa > sb) return 1;
    if (sa < sb) return -1;
    
    // Both zero
    if (sa == 0 && sb == 0) return 0;
    
    // Same sign - need to compare exponents and mantissas
    // Use subtraction: sign of (a - b) gives the answer
    int64_t s0[W + 10];
    s0[detail::IDX_ALLOCATED] = W + 6;
    
    detail::sub<W>(a.data(), b.data(), s0, W);
    
    if (s0[detail::IDX_SIGN_LENGTH] < 0) {
        return -1;
    } else if (s0[detail::IDX_SIGN_LENGTH] == 0) {
        return 0;
    } else {
        return 1;
    }
}

/**
 * @brief Fast comparison using raw data arrays
 * 
 * Internal function for use by detail namespace operations.
 * 
 * @tparam WORDS Precision in mantissa words
 * @param a First operand data array
 * @param b Second operand data array
 * @param mpnw Working precision
 * @return -1, 0, or +1 indicating comparison result
 */
namespace detail {

template <int WORDS>
KOKKOS_INLINE_FUNCTION
int cmp(const int64_t* a, const int64_t* b, int mpnw) {
    int sa = extract_sign(a[IDX_SIGN_LENGTH]);
    int sb = extract_sign(b[IDX_SIGN_LENGTH]);
    
    // Quick checks based on sign
    if (sa > sb) return 1;
    if (sa < sb) return -1;
    if (sa == 0 && sb == 0) return 0;
    
    // Same sign, non-zero - subtract and check result sign
    int64_t s0[WORDS + 10];
    s0[IDX_ALLOCATED] = mpnw + 6;
    
    sub<WORDS>(a, b, s0, mpnw);
    
    if (s0[IDX_SIGN_LENGTH] < 0) {
        return -1;
    } else if (s0[IDX_SIGN_LENGTH] == 0) {
        return 0;
    } else {
        return 1;
    }
}

/**
 * @brief Fast equality check
 * 
 * Optimized check for equality that can short-circuit on exponents.
 * 
 * @tparam WORDS Precision in mantissa words  
 * @param a First operand data array
 * @param b Second operand data array
 * @param mpnw Working precision
 * @return true if equal, false otherwise
 */
template <int WORDS>
KOKKOS_INLINE_FUNCTION
bool eq(const int64_t* a, const int64_t* b, int mpnw) {
    int sa = extract_sign(a[IDX_SIGN_LENGTH]);
    int sb = extract_sign(b[IDX_SIGN_LENGTH]);
    int na = extract_length(a[IDX_SIGN_LENGTH]);
    int nb = extract_length(b[IDX_SIGN_LENGTH]);
    
    // Quick checks
    if (sa != sb) return false;
    if (sa == 0 && sb == 0) return true;
    
    // Different exponents means different values
    if (a[IDX_EXPONENT] != b[IDX_EXPONENT]) return false;
    
    // Compare mantissa words (only up to the shorter length for now)
    int nmin = (na < nb) ? na : nb;
    int nmax = (na > nb) ? na : nb;
    
    for (int i = 0; i < nmin; ++i) {
        if (a[IDX_MANTISSA_START + i] != b[IDX_MANTISSA_START + i]) return false;
    }
    
    // If different lengths, check if extra words are all zero
    if (na > nb) {
        for (int i = nb; i < na; ++i) {
            if (a[IDX_MANTISSA_START + i] != 0) return false;
        }
    } else if (nb > na) {
        for (int i = na; i < nb; ++i) {
            if (b[IDX_MANTISSA_START + i] != 0) return false;
        }
    }
    
    return true;
}

/**
 * @brief Less-than comparison
 * 
 * @tparam WORDS Precision in mantissa words
 * @param a First operand data array
 * @param b Second operand data array
 * @param mpnw Working precision
 * @return true if a < b
 */
template <int WORDS>
KOKKOS_INLINE_FUNCTION
bool lt(const int64_t* a, const int64_t* b, int mpnw) {
    return cmp<WORDS>(a, b, mpnw) < 0;
}

/**
 * @brief Less-than-or-equal comparison
 * 
 * @tparam WORDS Precision in mantissa words
 * @param a First operand data array
 * @param b Second operand data array
 * @param mpnw Working precision
 * @return true if a <= b
 */
template <int WORDS>
KOKKOS_INLINE_FUNCTION
bool le(const int64_t* a, const int64_t* b, int mpnw) {
    return cmp<WORDS>(a, b, mpnw) <= 0;
}

/**
 * @brief Greater-than comparison
 * 
 * @tparam WORDS Precision in mantissa words
 * @param a First operand data array
 * @param b Second operand data array
 * @param mpnw Working precision
 * @return true if a > b
 */
template <int WORDS>
KOKKOS_INLINE_FUNCTION
bool gt(const int64_t* a, const int64_t* b, int mpnw) {
    return cmp<WORDS>(a, b, mpnw) > 0;
}

/**
 * @brief Greater-than-or-equal comparison
 * 
 * @tparam WORDS Precision in mantissa words
 * @param a First operand data array
 * @param b Second operand data array
 * @param mpnw Working precision
 * @return true if a >= b
 */
template <int WORDS>
KOKKOS_INLINE_FUNCTION
bool ge(const int64_t* a, const int64_t* b, int mpnw) {
    return cmp<WORDS>(a, b, mpnw) >= 0;
}

} // namespace detail
} // namespace mpfun

#endif // MPFUN_CORE_COMPARE_HPP
