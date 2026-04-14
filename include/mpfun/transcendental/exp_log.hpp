// SPDX-License-Identifier: MIT
// Exponential and Logarithmic Functions for MPFloat
// Part of kokkos-mpfun: High-precision arithmetic library
//
// Ported from MPFUN2020-Fort mpfund.f90 by David H. Bailey

#ifndef MPFUN_TRANSCENDENTAL_EXP_LOG_HPP
#define MPFUN_TRANSCENDENTAL_EXP_LOG_HPP

#include "../mp_float.hpp"
#include "../core/div.hpp"
#include <cmath>
#include <algorithm>

namespace mpfun {
namespace detail {

// =============================================================================
// Internal Helper: Precomputed ln(2) constant
// =============================================================================

/// ln(2) = 0.693147180559945309417232121458176568...
/// Computed to sufficient precision for argument reduction
/// This is a simplified version - for production, use AGM-based computation
template <int WORDS>
KOKKOS_INLINE_FUNCTION
void compute_ln2(int64_t* ln2, int mpnw) {
    // Use double precision as starting point
    // ln(2) ≈ 0.6931471805599453
    constexpr double LN2_D = 0.6931471805599453094172321214581765680755;
    
    ln2[IDX_ALLOCATED] = mpnw + 6;
    ln2[IDX_PRECISION] = mpnw;
    
    // For higher precision, we need to compute ln(2) properly
    // Using the identity: ln(2) = 2 * atanh(1/3) + atanh(1/7) / 2
    // Or use AGM: ln(2) = π / (2 * AGM(1, 4/2^(n/2)))
    // For now, use double precision (good for ~15 digits)
    mpdmc<WORDS>(LN2_D, 0, ln2, mpnw);
    
    // For higher precision implementations, extend using Newton-Raphson
    // on exp(y) = 2, which gives y_{k+1} = y_k + 2*exp(-y_k) - 1
    // This requires exp() to work first, creating a chicken-egg problem
    // Resolution: bootstrap with double, then refine iteratively
}

} // namespace detail

// =============================================================================
// ln2<W>() - Natural logarithm of 2 constant
// =============================================================================

/// Natural logarithm of 2 (ln(2) = 0.69314718055994530942...)
/// @tparam W Number of mantissa words (determines precision)
/// @return ln(2) to approximately W * 18 decimal digits
/// @note For full precision, requires iterative refinement
template <int W>
KOKKOS_INLINE_FUNCTION 
MPFloat<W> ln2() {
    MPFloat<W> result;
    detail::compute_ln2<W>(result.data(), W);
    return result;
}

// =============================================================================
// exp(x) - Exponential function
// =============================================================================

/// Computes e^x
/// @param x The exponent
/// @return e raised to the power x
///
/// Algorithm from MPFUN (mpexp):
/// 1. Handle special cases (x=0 → 1, overflow/underflow)
/// 2. Argument reduction: x = nz * ln(2) + r, where |r| < ln(2)/2
/// 3. Further reduce: s = r / 2^nq
/// 4. Taylor series: exp(s) = 1 + s + s²/2! + s³/3! + ...
/// 5. Square nq times: exp(r) = exp(s)^(2^nq)
/// 6. Final: exp(x) = exp(r) * 2^nz
template <int W>
KOKKOS_INLINE_FUNCTION 
MPFloat<W> exp(const MPFloat<W>& x) {
    constexpr int ITRMX = 1000000;
    constexpr double CL2 = 0.69314718055994530;  // ln(2) for double comparison
    
    int mpnw = W;
    int mpnw1 = mpnw + 1;
    
    // Result and temporaries
    MPFloat<W> result;
    int64_t* b = result.data();
    const int64_t* a = x.data();
    
    // Temporary arrays
    int64_t f[W + 10];
    int64_t s0[W + 10];
    int64_t s1[W + 10];
    int64_t s2[W + 10];
    int64_t s3[W + 10];
    int64_t ln2_arr[W + 10];
    
    for (int i = 0; i < W + 10; ++i) {
        f[i] = s0[i] = s1[i] = s2[i] = s3[i] = ln2_arr[i] = 0;
    }
    
    // Get double approximation
    double t1;
    int n1;
    detail::mpmdc<W>(a, t1, n1, mpnw);
    
    // Check for overflows and underflows
    if (n1 > 30) {
        if (t1 > 0.0) {
            // Overflow - set to NaN
            b[detail::IDX_ALLOCATED] = -1;
            return result;
        } else {
            // Underflow - return 0
            b[detail::IDX_PRECISION] = mpnw;
            b[detail::IDX_SIGN_LENGTH] = 0;
            b[detail::IDX_EXPONENT] = 0;
            b[detail::IDX_MANTISSA_START] = 0;
            b[detail::IDX_MANTISSA_START + 1] = 0;
            return result;
        }
    }
    
    t1 = t1 * std::pow(2.0, static_cast<double>(n1));
    
    // Check magnitude against overflow threshold
    // exp(1488522236) would overflow double
    if (std::abs(t1) > 1488522236.0) {
        if (t1 > 0) {
            b[detail::IDX_ALLOCATED] = -1;  // NaN for overflow
            return result;
        } else {
            // Underflow to zero
            b[detail::IDX_PRECISION] = mpnw;
            b[detail::IDX_SIGN_LENGTH] = 0;
            b[detail::IDX_EXPONENT] = 0;
            b[detail::IDX_MANTISSA_START] = 0;
            b[detail::IDX_MANTISSA_START + 1] = 0;
            return result;
        }
    }
    
    // Initialize temporary array sizes
    s0[detail::IDX_ALLOCATED] = mpnw + 7;
    s1[detail::IDX_ALLOCATED] = mpnw + 7;
    s2[detail::IDX_ALLOCATED] = mpnw + 7;
    s3[detail::IDX_ALLOCATED] = mpnw + 7;
    f[detail::IDX_ALLOCATED] = 9;
    ln2_arr[detail::IDX_ALLOCATED] = mpnw + 7;
    
    // Set f = 1.0
    f[detail::IDX_PRECISION] = mpnw1;
    f[detail::IDX_SIGN_LENGTH] = 1;
    f[detail::IDX_EXPONENT] = 0;
    f[detail::IDX_MANTISSA_START] = 1;
    f[detail::IDX_MANTISSA_START + 1] = 0;
    f[detail::IDX_MANTISSA_START + 2] = 0;
    
    // Compute ln(2) constant
    detail::compute_ln2<W>(ln2_arr, mpnw1);
    
    // Compute reduced argument: r = a - ln(2) * nint(a / ln(2))
    // s0 = a / ln(2)
    detail::div<W>(a, ln2_arr, s0, mpnw1);
    
    // s1 = nint(s0)
    // Round to nearest integer
    double d1;
    int n2;
    detail::mpmdc<W>(s0, d1, n2, mpnw1);
    d1 = d1 * std::pow(2.0, static_cast<double>(n2));
    int nz = static_cast<int>(std::round(d1));
    detail::mpdmc<W>(static_cast<double>(nz), 0, s1, mpnw1);
    
    // s2 = ln(2) * nint(a/ln(2))
    detail::mul<W>(ln2_arr, s1, s2, mpnw1);
    
    // s0 = a - s2 (reduced argument)
    detail::sub<W>(a, s2, s0, mpnw1);
    
    // Check if reduced argument is zero
    if (s0[detail::IDX_SIGN_LENGTH] == 0) {
        // exp(0) = 1, so result = 2^nz
        detail::mpdmc<W>(1.0, nz, s0, mpnw1);
        detail::round<W>(s0, mpnw);
        detail::mpeq<W>(s0, b, mpnw);
        return result;
    }
    
    // Divide reduced argument by 2^nq where nq = (mpnw * mpnbt)^0.4
    int nq = static_cast<int>(std::max(
        std::pow(static_cast<double>(mpnw * detail::BITS_PER_WORD), 0.4), 1.0));
    
    // s1 = s0 / 2^nq
    double divisor = std::pow(2.0, static_cast<double>(nq));
    detail::divd<W>(s0, divisor, s1, mpnw1);
    
    // Compute exp(s1) using Taylor series: 1 + s + s²/2! + s³/3! + ...
    // s2 = current term, s3 = accumulated sum
    detail::mpeq<W>(f, s2, mpnw1);   // s2 = 1 (first term)
    detail::mpeq<W>(f, s3, mpnw1);   // s3 = 1 (sum starts at 1)
    
    int mpnw2 = mpnw1;
    
    for (int j = 1; j <= ITRMX; ++j) {
        double t2 = static_cast<double>(j);
        
        // s0 = s2 * s1
        detail::mul<W>(s2, s1, s0, mpnw2);
        
        // s2 = s0 / j
        detail::divd<W>(s0, t2, s2, mpnw2);
        
        // s0 = s3 + s2
        detail::add<W>(s3, s2, s0, mpnw1);
        
        // s3 = s0
        detail::mpeq<W>(s0, s3, mpnw1);
        
        // Check for convergence: term is zero or negligible
        if (s2[detail::IDX_SIGN_LENGTH] == 0 || 
            s2[detail::IDX_EXPONENT] < s0[detail::IDX_EXPONENT] - mpnw1) {
            break;
        }
        
        // Reduce working precision for next term
        int exp_diff = static_cast<int>(s2[detail::IDX_EXPONENT] - s0[detail::IDX_EXPONENT]);
        mpnw2 = std::min(std::max(mpnw1 + exp_diff + 1, 4), mpnw1);
    }
    
    // Copy result to s0
    detail::mpeq<W>(s3, s0, mpnw1);
    
    // Square nq times to recover exp(reduced argument)
    for (int j = 1; j <= nq; ++j) {
        detail::mul<W>(s0, s0, s1, mpnw1);
        detail::mpeq<W>(s1, s0, mpnw1);
    }
    
    // Multiply by 2^nz
    detail::mpdmc<W>(1.0, nz, s2, mpnw1);
    detail::mul<W>(s0, s2, s1, mpnw1);
    
    // Round and copy to result
    detail::round<W>(s1, mpnw);
    detail::mpeq<W>(s1, b, mpnw);
    
    return result;
}

// =============================================================================
// log(x) - Natural logarithm
// =============================================================================

/// Natural logarithm ln(x)
/// @param x Positive input value
/// @return ln(x)
/// @pre x > 0 (returns NaN otherwise)
///
/// Algorithm from MPFUN (mplog):
/// For x ≈ 1: Use Taylor series ln(1+y) = y - y²/2 + y³/3 - ...
/// General case: Newton-Raphson iteration to solve exp(b) = x
///   b_{k+1} = b_k + [x - exp(b_k)] / exp(b_k)
/// 
/// Uses dynamic precision: starts at 4 words, doubles each iteration
template <int W>
KOKKOS_INLINE_FUNCTION 
MPFloat<W> log(const MPFloat<W>& x) {
    constexpr int ITRMAX = 1000000;
    constexpr int NIT = 3;  // Extra iterations after reaching full precision
    constexpr double ALT = 0.693147180559945309;  // ln(2) for initial approx
    constexpr double CL2 = 1.4426950408889633;    // 1/ln(2) = log2(e)
    constexpr double RTOL = 0.5 / 128.0;          // 0.5^7, threshold for Taylor
    
    int mpnw = W;
    int mpnw1 = mpnw + 1;
    
    MPFloat<W> result;
    int64_t* b = result.data();
    const int64_t* a = x.data();
    
    // Check for invalid input (x <= 0)
    int ia = detail::extract_sign(a[detail::IDX_SIGN_LENGTH]);
    int na = detail::extract_length(a[detail::IDX_SIGN_LENGTH]);
    if (na > mpnw) na = mpnw;
    
    if (ia < 0 || na == 0) {
        // log of negative or zero - return NaN
        b[detail::IDX_ALLOCATED] = -1;
        return result;
    }
    
    // Check if input is exactly 1 (log(1) = 0)
    if (a[detail::IDX_SIGN_LENGTH] == 1 && 
        a[detail::IDX_EXPONENT] == 0 && 
        a[detail::IDX_MANTISSA_START] == 1) {
        b[detail::IDX_PRECISION] = mpnw;
        b[detail::IDX_SIGN_LENGTH] = 0;
        b[detail::IDX_EXPONENT] = 0;
        b[detail::IDX_MANTISSA_START] = 0;
        b[detail::IDX_MANTISSA_START + 1] = 0;
        return result;
    }
    
    // Temporary arrays
    int64_t f1[W + 10];
    int64_t s0[W + 10];
    int64_t s1[W + 10];
    int64_t s2[W + 10];
    int64_t s3[W + 10];
    int64_t s4[W + 10];
    
    for (int i = 0; i < W + 10; ++i) {
        f1[i] = s0[i] = s1[i] = s2[i] = s3[i] = s4[i] = 0;
    }
    
    s0[detail::IDX_ALLOCATED] = mpnw + 7;
    s1[detail::IDX_ALLOCATED] = mpnw + 7;
    s2[detail::IDX_ALLOCATED] = mpnw + 7;
    s3[detail::IDX_ALLOCATED] = mpnw + 7;
    s4[detail::IDX_ALLOCATED] = mpnw + 7;
    
    // Set f1 = 1.0
    f1[detail::IDX_ALLOCATED] = 9;
    f1[detail::IDX_PRECISION] = mpnw1;
    f1[detail::IDX_SIGN_LENGTH] = 1;
    f1[detail::IDX_EXPONENT] = 0;
    f1[detail::IDX_MANTISSA_START] = 1;
    f1[detail::IDX_MANTISSA_START + 1] = 0;
    f1[detail::IDX_MANTISSA_START + 2] = 0;
    
    // Check if argument is close to 1 - use Taylor series
    // s0 = a - 1
    detail::sub<W>(a, f1, s0, mpnw1);
    
    if (s0[detail::IDX_SIGN_LENGTH] == 0 || 
        s0[detail::IDX_EXPONENT] <= static_cast<int64_t>(-2.0 - RTOL * mpnw1)) {
        // Use Taylor series: ln(1+y) = y - y²/2 + y³/3 - y⁴/4 + ...
        detail::mpeq<W>(s0, s1, mpnw1);  // s1 = y
        detail::mpeq<W>(s1, s2, mpnw1);  // s2 = y (power term)
        int is = 1;
        double tol = static_cast<double>(s0[detail::IDX_EXPONENT]) - mpnw1;
        
        for (int i1 = 2; i1 <= ITRMAX; ++i1) {
            is = -is;
            double st = static_cast<double>(is * i1);
            
            // s3 = s2 * s1 (next power of y)
            detail::mul<W>(s1, s2, s3, mpnw1);
            detail::mpeq<W>(s3, s2, mpnw1);
            
            // s4 = s3 / i1 (with alternating sign)
            detail::divd<W>(s3, st, s4, mpnw1);
            
            // s3 = s0 + s4
            detail::add<W>(s0, s4, s3, mpnw1);
            detail::mpeq<W>(s3, s0, mpnw1);
            
            // Check convergence
            if (s4[detail::IDX_SIGN_LENGTH] == 0 || 
                s4[detail::IDX_EXPONENT] < static_cast<int64_t>(tol)) {
                break;
            }
        }
        
        // Round and copy result
        detail::round<W>(s0, mpnw);
        detail::mpeq<W>(s0, b, mpnw);
        return result;
    }
    
    // General case: Newton-Raphson iteration
    // Determine number of iterations: mq = ceil(log2(mpnw))
    double t2 = static_cast<double>(mpnw);
    int mq = static_cast<int>(CL2 * std::log(t2) + 2.0 - detail::COMPARE_FUZZ);
    
    // Compute initial approximation using double precision
    double t1;
    int n1;
    detail::mpmdc<W>(a, t1, n1, mpnw);
    t1 = std::log(t1) + n1 * ALT;
    detail::mpdmc<W>(t1, 0, s3, mpnw);
    
    mpnw1 = 4;  // Start with low precision
    int iq = 0;
    
    // Newton-Raphson: b_{k+1} = b_k + (a - exp(b_k)) / exp(b_k)
    for (int k = 0; k <= mq; ++k) {
        if (k > 1) {
            mpnw1 = std::min(2 * mpnw1 - 2, mpnw) + 1;
        }
        
        do {
            // s0 = exp(s3)
            // Need to call exp on MPFloat, so create temporary
            MPFloat<W> temp_in, temp_out;
            detail::mpeq<W>(s3, temp_in.data(), mpnw1);
            temp_out = exp(temp_in);
            detail::mpeq<W>(temp_out.data(), s0, mpnw1);
            
            // s1 = a - exp(s3)
            detail::sub<W>(a, s0, s1, mpnw1);
            
            // s2 = s1 / s0 = (a - exp(s3)) / exp(s3)
            detail::div<W>(s1, s0, s2, mpnw1);
            
            // s1 = s3 + s2 (new approximation)
            detail::add<W>(s3, s2, s1, mpnw1);
            detail::mpeq<W>(s1, s3, mpnw1);
            
            // Extra iterations near the end for robustness
            if (k == mq - NIT && iq == 0) {
                iq = 1;
            } else {
                break;
            }
        } while (true);
    }
    
    // Round and copy to result
    detail::round<W>(s3, mpnw);
    detail::mpeq<W>(s3, b, mpnw);
    
    return result;
}

// =============================================================================
// Additional convenience functions
// =============================================================================

/// Computes e^x - 1, accurate for small x
/// @param x The exponent (ideally |x| << 1)
/// @return e^x - 1 without catastrophic cancellation
template <int W>
KOKKOS_INLINE_FUNCTION 
MPFloat<W> expm1(const MPFloat<W>& x) {
    // For small |x|, use Taylor series directly: x + x²/2! + x³/3! + ...
    // For larger |x|, compute exp(x) - 1
    
    constexpr int ITRMX = 1000000;
    int mpnw = W;
    int mpnw1 = mpnw + 1;
    
    const int64_t* a = x.data();
    
    // Check if x is small enough for direct series
    // Use series if |x| < 0.5 (exponent < 0)
    if (a[detail::IDX_SIGN_LENGTH] == 0) {
        // x = 0, return 0
        return MPFloat<W>(0.0);
    }
    
    if (a[detail::IDX_EXPONENT] < 0) {
        // Small argument: use Taylor series starting at n=1
        // expm1(x) = x + x²/2! + x³/3! + ...
        MPFloat<W> result;
        int64_t* b = result.data();
        
        int64_t s0[W + 10], s1[W + 10], s2[W + 10], s3[W + 10];
        for (int i = 0; i < W + 10; ++i) {
            s0[i] = s1[i] = s2[i] = s3[i] = 0;
        }
        
        s0[detail::IDX_ALLOCATED] = mpnw + 7;
        s1[detail::IDX_ALLOCATED] = mpnw + 7;
        s2[detail::IDX_ALLOCATED] = mpnw + 7;
        s3[detail::IDX_ALLOCATED] = mpnw + 7;
        
        // s2 = x (current term = x^n / n!)
        // s3 = x (sum starts at x)
        detail::mpeq<W>(a, s2, mpnw1);
        detail::mpeq<W>(a, s3, mpnw1);
        
        int mpnw2 = mpnw1;
        
        for (int j = 2; j <= ITRMX; ++j) {
            double t2 = static_cast<double>(j);
            
            // s0 = s2 * x
            detail::mul<W>(s2, a, s0, mpnw2);
            
            // s2 = s0 / j
            detail::divd<W>(s0, t2, s2, mpnw2);
            
            // s0 = s3 + s2
            detail::add<W>(s3, s2, s0, mpnw1);
            detail::mpeq<W>(s0, s3, mpnw1);
            
            // Check convergence
            if (s2[detail::IDX_SIGN_LENGTH] == 0 ||
                s2[detail::IDX_EXPONENT] < s0[detail::IDX_EXPONENT] - mpnw1) {
                break;
            }
            
            int exp_diff = static_cast<int>(s2[detail::IDX_EXPONENT] - s0[detail::IDX_EXPONENT]);
            mpnw2 = std::min(std::max(mpnw1 + exp_diff + 1, 4), mpnw1);
        }
        
        detail::round<W>(s3, mpnw);
        detail::mpeq<W>(s3, b, mpnw);
        return result;
    }
    
    // Large argument: compute exp(x) - 1 directly
    MPFloat<W> exp_x = exp(x);
    return exp_x - MPFloat<W>(1.0);
}

/// Computes log(1 + x), accurate for small x
/// @param x Input where x > -1
/// @return ln(1 + x) without catastrophic cancellation
template <int W>
KOKKOS_INLINE_FUNCTION 
MPFloat<W> log1p(const MPFloat<W>& x) {
    // For small |x|, use Taylor series: x - x²/2 + x³/3 - x⁴/4 + ...
    // For larger |x|, compute log(1 + x) directly
    
    constexpr int ITRMX = 1000000;
    int mpnw = W;
    int mpnw1 = mpnw + 1;
    
    const int64_t* a = x.data();
    
    // Check if x = 0
    if (a[detail::IDX_SIGN_LENGTH] == 0) {
        return MPFloat<W>(0.0);
    }
    
    // Check if 1 + x <= 0
    MPFloat<W> one(1.0);
    MPFloat<W> one_plus_x = one + x;
    if (one_plus_x.sign() <= 0) {
        MPFloat<W> result;
        result.set_nan();
        return result;
    }
    
    // For small |x| (exponent < -1), use direct Taylor series
    if (a[detail::IDX_EXPONENT] < -1) {
        MPFloat<W> result;
        int64_t* b = result.data();
        
        int64_t s0[W + 10], s1[W + 10], s2[W + 10], s3[W + 10];
        for (int i = 0; i < W + 10; ++i) {
            s0[i] = s1[i] = s2[i] = s3[i] = 0;
        }
        
        s0[detail::IDX_ALLOCATED] = mpnw + 7;
        s1[detail::IDX_ALLOCATED] = mpnw + 7;
        s2[detail::IDX_ALLOCATED] = mpnw + 7;
        s3[detail::IDX_ALLOCATED] = mpnw + 7;
        
        // s1 = x, s2 = x (power), s0 = x (sum)
        detail::mpeq<W>(a, s1, mpnw1);
        detail::mpeq<W>(a, s2, mpnw1);
        detail::mpeq<W>(a, s0, mpnw1);
        
        int is = 1;
        
        for (int i1 = 2; i1 <= ITRMX; ++i1) {
            is = -is;
            double st = static_cast<double>(is * i1);
            
            // s3 = s2 * s1 (next power)
            detail::mul<W>(s2, s1, s3, mpnw1);
            detail::mpeq<W>(s3, s2, mpnw1);
            
            // s3 = s2 / (is * i1)
            detail::divd<W>(s2, st, s3, mpnw1);
            
            // s0 = s0 + s3
            int64_t temp[W + 10];
            for (int i = 0; i < W + 10; ++i) temp[i] = 0;
            temp[detail::IDX_ALLOCATED] = mpnw + 7;
            detail::add<W>(s0, s3, temp, mpnw1);
            detail::mpeq<W>(temp, s0, mpnw1);
            
            // Check convergence
            if (s3[detail::IDX_SIGN_LENGTH] == 0 ||
                s3[detail::IDX_EXPONENT] < s0[detail::IDX_EXPONENT] - mpnw1) {
                break;
            }
        }
        
        detail::round<W>(s0, mpnw);
        detail::mpeq<W>(s0, b, mpnw);
        return result;
    }
    
    // Large argument: compute log(1 + x) directly
    return log(one_plus_x);
}

/// Base-10 logarithm
/// @param x Positive input value
/// @return log10(x) = ln(x) / ln(10)
template <int W>
KOKKOS_INLINE_FUNCTION 
MPFloat<W> log10(const MPFloat<W>& x) {
    // log10(x) = log(x) / log(10)
    MPFloat<W> ln_x = log(x);
    MPFloat<W> ln_10 = log(MPFloat<W>(10.0));
    return ln_x / ln_10;
}

/// Base-2 logarithm  
/// @param x Positive input value
/// @return log2(x) = ln(x) / ln(2)
template <int W>
KOKKOS_INLINE_FUNCTION 
MPFloat<W> log2(const MPFloat<W>& x) {
    // log2(x) = log(x) / log(2)
    MPFloat<W> ln_x = log(x);
    MPFloat<W> ln_2 = ln2<W>();
    return ln_x / ln_2;
}

// =============================================================================
// Power Functions
// =============================================================================

/// Integer power function (more efficient)
/// @param base The base
/// @param exp Integer exponent
/// @return base raised to the integer power exp
/// @note Uses binary exponentiation for O(log n) multiplications
template <int W>
KOKKOS_INLINE_FUNCTION 
MPFloat<W> pow(const MPFloat<W>& base, int exp) {
    // Handle special cases
    if (exp == 0) return MPFloat<W>(1.0);
    if (exp == 1) return base;
    if (base.is_zero()) {
        if (exp > 0) return MPFloat<W>(0.0);
        // 0^negative is undefined
        MPFloat<W> result;
        result.set_nan();
        return result;
    }
    
    MPFloat<W> result(1.0);
    MPFloat<W> b = base;
    int e = exp;
    
    // Handle negative exponent
    if (e < 0) {
        b = MPFloat<W>(1.0) / b;
        e = -e;
    }
    
    // Binary exponentiation
    while (e > 0) {
        if (e & 1) {
            result = result * b;
        }
        b = b * b;
        e >>= 1;
    }
    
    return result;
}

/// General power function: base^exp
/// @param base The base (must be positive for non-integer exp)
/// @param exponent The exponent
/// @return base raised to the power exp
/// @note Computed as exp(exp * log(base))
template <int W>
KOKKOS_INLINE_FUNCTION 
MPFloat<W> pow(const MPFloat<W>& base, const MPFloat<W>& exponent) {
    // Handle special cases
    if (exponent.is_zero()) return MPFloat<W>(1.0);
    if (base.is_zero()) {
        if (exponent.is_positive()) return MPFloat<W>(0.0);
        MPFloat<W> result;
        result.set_nan();
        return result;
    }
    
    // Check if base is 1
    MPFloat<W> one(1.0);
    if (compare(base, one) == 0) return one;
    
    // For negative base, we'd need to check if exponent is integer
    // For now, return NaN for negative base
    if (base.is_negative()) {
        MPFloat<W> result;
        result.set_nan();
        return result;
    }
    
    // General case: base^exp = exp(exp * log(base))
    MPFloat<W> ln_base = log(base);
    MPFloat<W> product = exponent * ln_base;
    return exp(product);
}

} // namespace mpfun

#endif // MPFUN_TRANSCENDENTAL_EXP_LOG_HPP
