// SPDX-License-Identifier: MIT
// Hyperbolic Functions for MPFloat
// Part of kokkos-mpfun: High-precision arbitrary precision library
//
// Ported from MPFUN2020-Fort subroutine mpcsshr (David H. Bailey)
//
// Algorithm notes (from TRANSCENDENTAL_ALGORITHMS.md):
// - For small |x| (exponent < -1): Taylor series for sinh, then cosh = sqrt(1 + sinh^2)
// - For larger |x|: exponential definition sinh = (e^x - e^-x)/2, cosh = (e^x + e^-x)/2

#ifndef MPFUN_TRANSCENDENTAL_HYPERBOLIC_HPP
#define MPFUN_TRANSCENDENTAL_HYPERBOLIC_HPP

#include "../core/representation.hpp"
#include "../core/add.hpp"
#include "../core/mul.hpp"
#include "../core/div.hpp"
#include "../core/sqrt.hpp"

// Forward declaration of MPFloat to avoid circular includes
namespace mpfun {
template <int WORDS> class MPFloat;

template <int W>
KOKKOS_INLINE_FUNCTION MPFloat<W> abs(const MPFloat<W>& x);
}

namespace mpfun {

// Forward declarations
template <int W>
KOKKOS_INLINE_FUNCTION 
void sinhcosh(const MPFloat<W>& x, MPFloat<W>& sinh_out, MPFloat<W>& cosh_out);

template <int W>
KOKKOS_INLINE_FUNCTION 
MPFloat<W> sinh(const MPFloat<W>& x);

template <int W>
KOKKOS_INLINE_FUNCTION 
MPFloat<W> cosh(const MPFloat<W>& x);

namespace detail {

/**
 * @brief Compute hyperbolic sine and cosine together (internal implementation)
 * 
 * This implements the algorithm from mpcsshr in mpfund.f90:
 * 
 * For small arguments (exponent < -1):
 *   - Use Taylor series for sinh: x + x³/3! + x⁵/5! + ...
 *   - Compute cosh from sinh: cosh = sqrt(1 + sinh²)
 *   This avoids precision loss from exponential subtraction.
 * 
 * For larger arguments:
 *   - Use exponential: t = exp(x)
 *   - cosh = (t + 1/t) / 2
 *   - sinh = (t - 1/t) / 2
 * 
 * @tparam WORDS Precision in mantissa words
 * @param a Input value
 * @param x Output cosh(a)
 * @param y Output sinh(a)
 * @param mpnw Working precision in words
 */
template <int WORDS>
KOKKOS_INLINE_FUNCTION
void sinhcosh_impl(const int64_t* a, int64_t* x, int64_t* y, int mpnw) {
    constexpr int ITRMX = 1000000;
    
    // Temporary arrays
    int64_t s0[WORDS + 10];
    int64_t s1[WORDS + 10];
    int64_t s2[WORDS + 10];
    int64_t s3[WORDS + 10];
    int64_t f[10];  // For constant 1
    
    // Initialize temporaries
    for (int i = 0; i < WORDS + 10; ++i) {
        s0[i] = s1[i] = s2[i] = s3[i] = 0;
    }
    
    s0[IDX_ALLOCATED] = mpnw + 7;
    s1[IDX_ALLOCATED] = mpnw + 7;
    s2[IDX_ALLOCATED] = mpnw + 7;
    s3[IDX_ALLOCATED] = mpnw + 7;
    
    int mpnw1 = mpnw + 1;
    
    // Initialize f = 1.0
    f[IDX_ALLOCATED] = 9;
    f[IDX_PRECISION] = mpnw;
    f[IDX_SIGN_LENGTH] = 1;  // sign=1, length=1
    f[IDX_EXPONENT] = 0;
    f[IDX_MANTISSA_START] = 1;
    f[IDX_MANTISSA_START + 1] = 0;
    f[6] = 0;
    f[7] = 0;
    
    // Get exponent to decide which algorithm to use
    int64_t a_exp = a[IDX_EXPONENT];
    
    // If argument is very small (exponent < -1), use Taylor series for sinh
    // This avoids accuracy loss from exp(x) - exp(-x) when x is small
    if (a_exp < -1) {
        // Taylor series: sinh(x) = x + x³/3! + x⁵/5! + ...
        // s0 will accumulate the sum, s1 is the current term
        mpeq<WORDS>(a, s0, mpnw1);  // s0 = x (first term)
        mpeq<WORDS>(a, s1, mpnw1);  // s1 = x (current term)
        mul<WORDS>(s0, s0, s2, mpnw1);  // s2 = x²
        
        int mpnw2 = mpnw1;
        
        // Compute x + x³/3! + x⁵/5! + ...
        for (int j = 1; j <= ITRMX; ++j) {
            // t2 = (2j)(2j+1) - the factorial denominator factor
            double t2 = (2.0 * j) * (2.0 * j + 1.0);
            
            // s1 = s1 * x² / t2 (next term)
            mul<WORDS>(s2, s1, s3, mpnw2);       // s3 = s1 * x²
            divd<WORDS>(s3, t2, s1, mpnw2);      // s1 = s3 / t2
            
            // s0 = s0 + s1 (accumulate)
            add<WORDS>(s1, s0, s3, mpnw1);
            mpeq<WORDS>(s3, s0, mpnw1);
            
            // Check for convergence: term is zero or negligible compared to sum
            if (s1[IDX_SIGN_LENGTH] == 0 || 
                s1[IDX_EXPONENT] < s0[IDX_EXPONENT] - mpnw1) {
                break;
            }
            
            // Reduce working precision for next term (linear precision reduction)
            int prec_diff = static_cast<int>(s1[IDX_EXPONENT] - s0[IDX_EXPONENT]);
            mpnw2 = mpnw1 + prec_diff + 1;
            if (mpnw2 < 4) mpnw2 = 4;
            if (mpnw2 > mpnw1) mpnw2 = mpnw1;
        }
        
        // Now s0 = sinh(x)
        // Compute cosh(x) = sqrt(1 + sinh²(x))
        mul<WORDS>(s0, s0, s2, mpnw1);      // s2 = sinh²
        add<WORDS>(f, s2, s3, mpnw1);       // s3 = 1 + sinh²
        sqrt<WORDS>(s3, s1, mpnw1);         // s1 = sqrt(1 + sinh²) = cosh
        
        // Round and output
        round<WORDS>(s1, mpnw);
        mpeq<WORDS>(s1, x, mpnw);  // x = cosh
        round<WORDS>(s0, mpnw);
        mpeq<WORDS>(s0, y, mpnw);  // y = sinh
    } else {
        // For larger arguments, use exponential definition
        // t = exp(a)
        // cosh = (t + 1/t) / 2
        // sinh = (t - 1/t) / 2
        
        // s0 = exp(a)
        exp_impl<WORDS>(a, s0, mpnw1);
        
        // s1 = 1 / exp(a) = exp(-a)
        div<WORDS>(f, s0, s1, mpnw1);
        
        // s2 = exp(a) + exp(-a)
        add<WORDS>(s0, s1, s2, mpnw1);
        
        // s3 = (exp(a) + exp(-a)) / 2 = cosh
        muld<WORDS>(s2, 0.5, s3, mpnw1);
        round<WORDS>(s3, mpnw);
        mpeq<WORDS>(s3, x, mpnw);  // x = cosh
        
        // s2 = exp(a) - exp(-a)
        sub<WORDS>(s0, s1, s2, mpnw1);
        
        // s3 = (exp(a) - exp(-a)) / 2 = sinh
        muld<WORDS>(s2, 0.5, s3, mpnw1);
        round<WORDS>(s3, mpnw);
        mpeq<WORDS>(s3, y, mpnw);  // y = sinh
    }
}

/**
 * @brief Compute exponential using argument reduction + Taylor series
 * 
 * Uses the algorithm from mpexp in mpfund.f90:
 * 1. Argument reduction: x = k*log(2) + r where |r| < log(2)/2
 * 2. Further reduction: s = r / 2^NQ
 * 3. Taylor series for exp(s)
 * 4. Square NQ times to get exp(r)
 * 5. Multiply by 2^k
 * 
 * @tparam WORDS Precision in mantissa words
 * @param a Input value
 * @param b Output exp(a)
 * @param mpnw Working precision in words
 */
template <int WORDS>
KOKKOS_INLINE_FUNCTION
void exp_impl(const int64_t* a, int64_t* b, int mpnw) {
    constexpr int ITRMX = 1000000;
    constexpr double LOG2 = 0.693147180559945309417232121458176568;  // ln(2)
    
    int mpnw1 = mpnw + 1;
    
    // Temporary arrays
    int64_t s0[WORDS + 10];
    int64_t s1[WORDS + 10];
    int64_t s2[WORDS + 10];
    int64_t s3[WORDS + 10];
    int64_t f[10];
    
    for (int i = 0; i < WORDS + 10; ++i) {
        s0[i] = s1[i] = s2[i] = s3[i] = 0;
    }
    
    s0[IDX_ALLOCATED] = mpnw + 7;
    s1[IDX_ALLOCATED] = mpnw + 7;
    s2[IDX_ALLOCATED] = mpnw + 7;
    s3[IDX_ALLOCATED] = mpnw + 7;
    
    // Initialize f = 1.0
    f[IDX_ALLOCATED] = 9;
    f[IDX_PRECISION] = mpnw;
    f[IDX_SIGN_LENGTH] = 1;
    f[IDX_EXPONENT] = 0;
    f[IDX_MANTISSA_START] = 1;
    f[IDX_MANTISSA_START + 1] = 0;
    f[6] = 0;
    f[7] = 0;
    
    // Check for zero input: exp(0) = 1
    if (a[IDX_SIGN_LENGTH] == 0) {
        mpeq<WORDS>(f, b, mpnw);
        return;
    }
    
    // Get double approximation for argument reduction
    double da;
    int n_exp;
    mpmdc<WORDS>(a, da, n_exp, mpnw);
    da = da * Kokkos::pow(2.0, static_cast<double>(n_exp));
    
    // Argument reduction: a = nz * log(2) + t, where |t| < log(2)/2
    int nz = static_cast<int>(da / LOG2 + 0.5);
    if (da < 0.0) nz = static_cast<int>(da / LOG2 - 0.5);
    
    // t = a - nz * log(2)
    mpdmc<WORDS>(static_cast<double>(nz) * LOG2, 0, s1, mpnw1);
    sub<WORDS>(a, s1, s0, mpnw1);  // s0 = t = a - nz*log(2)
    
    // Further reduction: divide by 2^NQ
    // NQ ≈ (mpnw * bits_per_word)^0.4
    int nq = static_cast<int>(Kokkos::pow(static_cast<double>(mpnw) * 60.0, 0.4));
    if (nq < 1) nq = 1;
    
    // s1 = t / 2^NQ
    double t_scale = Kokkos::pow(2.0, static_cast<double>(-nq));
    muld<WORDS>(s0, t_scale, s1, mpnw1);
    
    // Taylor series: exp(s) = 1 + s + s²/2! + s³/3! + ...
    // s2 = s (current term), s0 = 1 (sum)
    mpeq<WORDS>(s1, s2, mpnw1);  // s2 = s (first term after 1)
    mpeq<WORDS>(f, s0, mpnw1);   // s0 = 1
    add<WORDS>(s0, s2, s3, mpnw1);  // s3 = 1 + s
    mpeq<WORDS>(s3, s0, mpnw1);     // s0 = 1 + s
    
    int mpnw2 = mpnw1;
    
    for (int j = 2; j <= ITRMX; ++j) {
        // s2 = s2 * s / j (next term)
        mul<WORDS>(s2, s1, s3, mpnw2);
        divd<WORDS>(s3, static_cast<double>(j), s2, mpnw2);
        
        // s0 = s0 + s2
        add<WORDS>(s0, s2, s3, mpnw1);
        mpeq<WORDS>(s3, s0, mpnw1);
        
        // Check convergence
        if (s2[IDX_SIGN_LENGTH] == 0 ||
            s2[IDX_EXPONENT] < s0[IDX_EXPONENT] - mpnw1) {
            break;
        }
        
        // Linear precision reduction
        int prec_diff = static_cast<int>(s2[IDX_EXPONENT] - s0[IDX_EXPONENT]);
        mpnw2 = mpnw1 + prec_diff + 1;
        if (mpnw2 < 4) mpnw2 = 4;
        if (mpnw2 > mpnw1) mpnw2 = mpnw1;
    }
    
    // Recovery: square NQ times to get exp(t) from exp(t/2^NQ)
    for (int j = 0; j < nq; ++j) {
        mul<WORDS>(s0, s0, s1, mpnw1);
        mpeq<WORDS>(s1, s0, mpnw1);
    }
    
    // Final adjustment: multiply by 2^nz
    // This is done by adjusting the exponent
    if (nz != 0) {
        s0[IDX_EXPONENT] = s0[IDX_EXPONENT] + nz;
    }
    
    round<WORDS>(s0, mpnw);
    mpeq<WORDS>(s0, b, mpnw);
}

/**
 * @brief Compute natural logarithm using Newton-Raphson iteration
 * 
 * Uses the algorithm from mplog in mpfund.f90:
 * For general case: Newton iteration to solve exp(b) = a
 *   x_{k+1} = x_k + [a - exp(x_k)] / exp(x_k)
 * 
 * For input close to 1: Taylor series for log(1+x)
 * 
 * @tparam WORDS Precision in mantissa words
 * @param a Input value (must be positive)
 * @param b Output log(a)
 * @param mpnw Working precision in words
 */
template <int WORDS>
KOKKOS_INLINE_FUNCTION
void log_impl(const int64_t* a, int64_t* b, int mpnw) {
    constexpr int NIT = 3;  // Extra iterations at full precision
    constexpr double LOG2 = 0.693147180559945309417232121458176568;
    
    int ia = extract_sign(a[IDX_SIGN_LENGTH]);
    int na = extract_length(a[IDX_SIGN_LENGTH]);
    
    // Check for invalid input (a <= 0)
    if (ia <= 0 || na == 0) {
        b[IDX_ALLOCATED] = -1;  // NaN
        b[IDX_PRECISION] = mpnw;
        b[IDX_SIGN_LENGTH] = 0;
        return;
    }
    
    int mpnw1 = mpnw + 1;
    
    // Temporary arrays
    int64_t s0[WORDS + 10];
    int64_t s1[WORDS + 10];
    int64_t s2[WORDS + 10];
    int64_t s3[WORDS + 10];
    int64_t s4[WORDS + 10];
    
    for (int i = 0; i < WORDS + 10; ++i) {
        s0[i] = s1[i] = s2[i] = s3[i] = s4[i] = 0;
    }
    
    s0[IDX_ALLOCATED] = mpnw + 7;
    s1[IDX_ALLOCATED] = mpnw + 7;
    s2[IDX_ALLOCATED] = mpnw + 7;
    s3[IDX_ALLOCATED] = mpnw + 7;
    s4[IDX_ALLOCATED] = mpnw + 7;
    
    // Get initial approximation using double precision
    double da;
    int n_exp;
    mpmdc<WORDS>(a, da, n_exp, mpnw);
    
    // Initial approximation: log(da * 2^n) = log(da) + n*log(2)
    double t1 = Kokkos::log(da) + static_cast<double>(n_exp) * LOG2;
    mpdmc<WORDS>(t1, 0, s2, mpnw1);
    
    // Newton-Raphson iteration with dynamic precision
    // x_{k+1} = x_k + [a - exp(x_k)] / exp(x_k)
    //         = x_k + a/exp(x_k) - 1
    //         = x_k + a*exp(-x_k) - 1
    
    // Number of iterations: mq = ceil(log2(mpnw))
    double t2 = static_cast<double>(mpnw);
    int mq = static_cast<int>(LOG2_E * Kokkos::log(t2) + 1.0 - COMPARE_FUZZ);
    
    int mpnw_iter = 5;
    int iq = 0;
    
    for (int k = 1; k <= mq; ++k) {
        if (k > 1) {
            mpnw_iter = 2 * mpnw_iter - 2;
            if (mpnw_iter > mpnw) mpnw_iter = mpnw;
            mpnw_iter += 1;
        }
        
        do {
            // s0 = exp(x_k)
            exp_impl<WORDS>(s2, s0, mpnw_iter);
            
            // s1 = a - exp(x_k)
            sub<WORDS>(a, s0, s1, mpnw_iter);
            
            // s3 = (a - exp(x_k)) / exp(x_k)
            div<WORDS>(s1, s0, s3, mpnw_iter);
            
            // s4 = x_k + correction
            add<WORDS>(s2, s3, s4, mpnw_iter);
            mpeq<WORDS>(s4, s2, mpnw_iter);
            
            // Extra iterations at full precision
            if (k == mq - NIT && iq == 0) {
                iq = 1;
            } else {
                break;
            }
        } while (true);
    }
    
    round<WORDS>(s2, mpnw);
    mpeq<WORDS>(s2, b, mpnw);
}

} // namespace detail

// =============================================================================
// Public API - Basic Hyperbolic Functions
// =============================================================================

/**
 * @brief Simultaneous hyperbolic sine and cosine (most efficient)
 * 
 * Computes both sinh(x) and cosh(x) together, sharing intermediate
 * computations. This is more efficient than calling sinh() and cosh()
 * separately when both values are needed.
 * 
 * @tparam W Precision in mantissa words
 * @param x Input value
 * @param sinh_out Output: sinh(x)
 * @param cosh_out Output: cosh(x)
 */
template <int W>
KOKKOS_INLINE_FUNCTION 
void sinhcosh(const MPFloat<W>& x, MPFloat<W>& sinh_out, MPFloat<W>& cosh_out) {
    detail::sinhcosh_impl<W>(x.data(), cosh_out.data(), sinh_out.data(), W);
}

/**
 * @brief Hyperbolic sine
 * 
 * Computes sinh(x) = (e^x - e^(-x)) / 2
 * 
 * For small |x|, uses Taylor series to avoid precision loss.
 * 
 * @tparam W Precision in mantissa words
 * @param x Input value
 * @return sinh(x)
 */
template <int W>
KOKKOS_INLINE_FUNCTION 
MPFloat<W> sinh(const MPFloat<W>& x) {
    MPFloat<W> sinh_out, cosh_out;
    sinhcosh(x, sinh_out, cosh_out);
    return sinh_out;
}

/**
 * @brief Hyperbolic cosine
 * 
 * Computes cosh(x) = (e^x + e^(-x)) / 2
 * 
 * Note: cosh(x) >= 1 for all real x, and cosh(-x) = cosh(x).
 * 
 * @tparam W Precision in mantissa words
 * @param x Input value
 * @return cosh(x)
 */
template <int W>
KOKKOS_INLINE_FUNCTION 
MPFloat<W> cosh(const MPFloat<W>& x) {
    MPFloat<W> sinh_out, cosh_out;
    sinhcosh(x, sinh_out, cosh_out);
    return cosh_out;
}

/**
 * @brief Hyperbolic tangent
 * 
 * Computes tanh(x) = sinh(x) / cosh(x)
 * 
 * Result is always in (-1, 1); approaches ±1 for large |x|.
 * 
 * @tparam W Precision in mantissa words
 * @param x Input value
 * @return tanh(x)
 */
template <int W>
KOKKOS_INLINE_FUNCTION 
MPFloat<W> tanh(const MPFloat<W>& x) {
    MPFloat<W> sinh_val, cosh_val;
    sinhcosh(x, sinh_val, cosh_val);
    return sinh_val / cosh_val;
}

// =============================================================================
// Public API - Inverse Hyperbolic Functions
// =============================================================================

/**
 * @brief Inverse hyperbolic sine
 * 
 * Computes asinh(x) = log(x + sqrt(x² + 1))
 * 
 * This function is defined for all real x and is an odd function.
 * 
 * @tparam W Precision in mantissa words
 * @param x Input value (any real number)
 * @return asinh(x)
 */
template <int W>
KOKKOS_INLINE_FUNCTION 
MPFloat<W> asinh(const MPFloat<W>& x) {
    // asinh(x) = log(x + sqrt(x² + 1))
    MPFloat<W> one(1);
    MPFloat<W> x_squared = x * x;
    MPFloat<W> sqrt_arg = x_squared + one;
    MPFloat<W> sqrt_val = sqrt(sqrt_arg);
    MPFloat<W> log_arg = x + sqrt_val;
    
    // Use log implementation
    MPFloat<W> result;
    detail::log_impl<W>(log_arg.data(), result.data(), W);
    return result;
}

/**
 * @brief Inverse hyperbolic cosine
 * 
 * Computes acosh(x) = log(x + sqrt(x² - 1)) for x >= 1
 * 
 * Returns NaN for x < 1.
 * 
 * @tparam W Precision in mantissa words
 * @param x Input value (must be >= 1)
 * @return acosh(x), or NaN if x < 1
 */
template <int W>
KOKKOS_INLINE_FUNCTION 
MPFloat<W> acosh(const MPFloat<W>& x) {
    // acosh(x) = log(x + sqrt(x² - 1)) for x >= 1
    MPFloat<W> one(1);
    
    // Check x >= 1
    if (x < one) {
        MPFloat<W> result;
        result.set_nan();
        return result;
    }
    
    MPFloat<W> x_squared = x * x;
    MPFloat<W> sqrt_arg = x_squared - one;
    MPFloat<W> sqrt_val = sqrt(sqrt_arg);
    MPFloat<W> log_arg = x + sqrt_val;
    
    MPFloat<W> result;
    detail::log_impl<W>(log_arg.data(), result.data(), W);
    return result;
}

/**
 * @brief Inverse hyperbolic tangent
 * 
 * Computes atanh(x) = (1/2) * log((1+x)/(1-x)) for |x| < 1
 * 
 * Returns ±inf at x = ±1, NaN for |x| > 1.
 * 
 * @tparam W Precision in mantissa words
 * @param x Input value (must satisfy |x| < 1)
 * @return atanh(x)
 */
template <int W>
KOKKOS_INLINE_FUNCTION 
MPFloat<W> atanh(const MPFloat<W>& x) {
    // atanh(x) = (1/2) * log((1+x)/(1-x)) for |x| < 1
    MPFloat<W> one(1);
    MPFloat<W> abs_x = abs(x);
    
    // Check |x| < 1
    if (abs_x >= one) {
        MPFloat<W> result;
        result.set_nan();  // Could also return ±inf at x = ±1
        return result;
    }
    
    // Compute (1+x)/(1-x)
    MPFloat<W> one_plus_x = one + x;
    MPFloat<W> one_minus_x = one - x;
    MPFloat<W> ratio = one_plus_x / one_minus_x;
    
    // log((1+x)/(1-x))
    MPFloat<W> log_val;
    detail::log_impl<W>(ratio.data(), log_val.data(), W);
    
    // Multiply by 0.5
    MPFloat<W> half(0.5);
    return log_val * half;
}

// =============================================================================
// Additional Utility Functions
// =============================================================================

/**
 * @brief Hyperbolic secant
 * 
 * Computes sech(x) = 1 / cosh(x)
 * 
 * @tparam W Precision in mantissa words
 * @param x Input value
 * @return sech(x)
 */
template <int W>
KOKKOS_INLINE_FUNCTION 
MPFloat<W> sech(const MPFloat<W>& x) {
    MPFloat<W> one(1);
    return one / cosh(x);
}

/**
 * @brief Hyperbolic cosecant
 * 
 * Computes csch(x) = 1 / sinh(x)
 * 
 * @tparam W Precision in mantissa words
 * @param x Input value (x != 0)
 * @return csch(x)
 */
template <int W>
KOKKOS_INLINE_FUNCTION 
MPFloat<W> csch(const MPFloat<W>& x) {
    MPFloat<W> one(1);
    return one / sinh(x);
}

/**
 * @brief Hyperbolic cotangent
 * 
 * Computes coth(x) = cosh(x) / sinh(x) = 1 / tanh(x)
 * 
 * @tparam W Precision in mantissa words
 * @param x Input value (x != 0)
 * @return coth(x)
 */
template <int W>
KOKKOS_INLINE_FUNCTION 
MPFloat<W> coth(const MPFloat<W>& x) {
    MPFloat<W> sinh_val, cosh_val;
    sinhcosh(x, sinh_val, cosh_val);
    return cosh_val / sinh_val;
}

} // namespace mpfun

#endif // MPFUN_TRANSCENDENTAL_HYPERBOLIC_HPP
