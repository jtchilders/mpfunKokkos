// SPDX-License-Identifier: MIT
// Exponential and Logarithmic Functions for MPFloat
// Part of kokkos-mpfun: High-precision arithmetic library

#ifndef MPFUN_TRANSCENDENTAL_EXP_LOG_HPP
#define MPFUN_TRANSCENDENTAL_EXP_LOG_HPP

#include "../mp_float.hpp"
#include <Kokkos_Core.hpp>

namespace mpfun {

// =============================================================================
// Exponential Functions
// =============================================================================

/// Computes e^x
/// @param x The exponent
/// @return e raised to the power x
/// @note Uses argument reduction: e^x = 2^k * e^r where |r| < ln(2)/2
///       Then Taylor series for e^r
template <int W>
KOKKOS_INLINE_FUNCTION 
MPFloat<W> exp(const MPFloat<W>& x) {
    // TODO: implement using algorithm from TRANSCENDENTAL_ALGORITHMS.md
    // Algorithm: Argument reduction + Taylor series
    // 1. Handle special cases (x=0 -> 1, x very negative -> 0)
    // 2. Reduce: x = k*ln(2) + r, |r| < ln(2)/2
    // 3. Compute e^r via Taylor: sum_{n=0}^{N} r^n/n!
    // 4. Return 2^k * e^r via exponent adjustment
    MPFloat<W> result;
    return result;
}

/// Computes e^x - 1, accurate for small x
/// @param x The exponent (ideally |x| << 1)
/// @return e^x - 1 without catastrophic cancellation
/// @note For |x| > threshold, falls back to exp(x) - 1
template <int W>
KOKKOS_INLINE_FUNCTION 
MPFloat<W> expm1(const MPFloat<W>& x) {
    // TODO: implement using algorithm from TRANSCENDENTAL_ALGORITHMS.md
    // Algorithm: Direct Taylor series for small |x|
    // For |x| < threshold (e.g., 0.5):
    //   sum_{n=1}^{N} x^n/n!  (note: starts at n=1)
    // For larger |x|: compute exp(x) - 1 directly
    MPFloat<W> result;
    return result;
}

// =============================================================================
// Logarithmic Functions
// =============================================================================

/// Natural logarithm ln(x)
/// @param x Positive input value
/// @return ln(x)
/// @pre x > 0 (returns NaN otherwise)
/// @note Uses argument reduction + AGM or Taylor series
template <int W>
KOKKOS_INLINE_FUNCTION 
MPFloat<W> log(const MPFloat<W>& x) {
    // TODO: implement using algorithm from TRANSCENDENTAL_ALGORITHMS.md
    // Algorithm options (choose based on precision needs):
    // 
    // Option A: AGM-based (quadratic convergence, better for high precision)
    //   log(x) = pi / (2 * AGM(1, 4/s)) - m*log(2)
    //   where s = x * 2^m chosen so s > 2^(p/2)
    //
    // Option B: Argument reduction + Taylor
    //   1. Write x = 2^k * f where 1 <= f < 2
    //   2. Let r = (f-1)/(f+1), so f = (1+r)/(1-r)
    //   3. log(f) = 2 * sum_{n=0}^{N} r^(2n+1)/(2n+1)
    //   4. Return k*ln(2) + log(f)
    MPFloat<W> result;
    return result;
}

/// Computes log(1 + x), accurate for small x
/// @param x Input where x > -1
/// @return ln(1 + x) without catastrophic cancellation
/// @pre x > -1 (returns NaN otherwise)
template <int W>
KOKKOS_INLINE_FUNCTION 
MPFloat<W> log1p(const MPFloat<W>& x) {
    // TODO: implement using algorithm from TRANSCENDENTAL_ALGORITHMS.md
    // Algorithm: For |x| < threshold:
    //   Use series: x - x^2/2 + x^3/3 - x^4/4 + ...
    //   Or: 2 * atanh(x/(2+x)) = 2 * atanh(r) where r = x/(2+x)
    // For larger |x|: compute log(1 + x) directly
    MPFloat<W> result;
    return result;
}

/// Base-10 logarithm
/// @param x Positive input value
/// @return log10(x) = ln(x) / ln(10)
/// @pre x > 0
template <int W>
KOKKOS_INLINE_FUNCTION 
MPFloat<W> log10(const MPFloat<W>& x) {
    // TODO: implement
    // Simply: log(x) / log(10)
    // Precompute 1/log(10) as a constant for efficiency
    MPFloat<W> result;
    return result;
}

/// Base-2 logarithm  
/// @param x Positive input value
/// @return log2(x) = ln(x) / ln(2)
/// @pre x > 0
template <int W>
KOKKOS_INLINE_FUNCTION 
MPFloat<W> log2(const MPFloat<W>& x) {
    // TODO: implement
    // Simply: log(x) / log(2)
    // Precompute 1/log(2) = log2(e) as a constant for efficiency
    // Can also extract exponent for fast integer part
    MPFloat<W> result;
    return result;
}

// =============================================================================
// Power Functions
// =============================================================================

/// General power function: base^exp
/// @param base The base (must be positive for non-integer exp)
/// @param exp The exponent
/// @return base raised to the power exp
/// @note Computed as exp(exp * log(base))
template <int W>
KOKKOS_INLINE_FUNCTION 
MPFloat<W> pow(const MPFloat<W>& base, const MPFloat<W>& exp) {
    // TODO: implement
    // Algorithm: base^exp = exp(exp * log(base))
    // Special cases:
    //   - base = 0: return 0 if exp > 0, NaN if exp <= 0
    //   - exp = 0: return 1
    //   - base < 0: NaN (unless exp is integer, then use integer version)
    //   - base = 1: return 1
    MPFloat<W> result;
    return result;
}

/// Integer power function (more efficient)
/// @param base The base
/// @param exp Integer exponent
/// @return base raised to the integer power exp
/// @note Uses binary exponentiation for O(log n) multiplications
template <int W>
KOKKOS_INLINE_FUNCTION 
MPFloat<W> pow(const MPFloat<W>& base, int exp) {
    // TODO: implement using binary exponentiation
    // Algorithm (exponentiation by squaring):
    //   result = 1
    //   if exp < 0: base = 1/base, exp = -exp
    //   while exp > 0:
    //     if exp is odd: result *= base
    //     base *= base
    //     exp >>= 1
    //   return result
    MPFloat<W> result;
    return result;
}

} // namespace mpfun

#endif // MPFUN_TRANSCENDENTAL_EXP_LOG_HPP
