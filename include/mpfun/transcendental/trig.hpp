// SPDX-License-Identifier: MIT
// Trigonometric Functions for MPFloat
// Part of kokkos-mpfun: High-precision arithmetic library
// 
// Implements sin, cos, sincos, tan, atan, atan2, asin, acos
// 
// Ported from MPFUN2020-Fort mpfund.f90 (mpcssn, mpang)

#ifndef MPFUN_TRANSCENDENTAL_TRIG_HPP
#define MPFUN_TRANSCENDENTAL_TRIG_HPP

#include "../core/representation.hpp"
#include "../core/add.hpp"
#include "../core/mul.hpp"
#include "../core/div.hpp"
#include "../core/sqrt.hpp"
#include <cmath>

namespace mpfun {

// Forward declarations
template <int W> class MPFloat;
template <int W> KOKKOS_INLINE_FUNCTION MPFloat<W> sqrt(const MPFloat<W>& x);

namespace detail {

// =============================================================================
// Helper: Round to nearest integer (mpnint)
// =============================================================================

/**
 * @brief Round MPFloat to nearest integer
 * 
 * Ported from mpnint in mpfunb.f90.
 * 
 * @tparam WORDS Precision in mantissa words
 * @param a Input MPFloat
 * @param b Output (rounded to nearest integer)
 * @param mpnw Working precision
 */
template <int WORDS>
KOKKOS_INLINE_FUNCTION
void mpnint(const int64_t* a, int64_t* b, int mpnw) {
    int ia = extract_sign(a[IDX_SIGN_LENGTH]);
    int na = extract_length(a[IDX_SIGN_LENGTH]);
    
    if (na == 0) {
        b[IDX_PRECISION] = mpnw;
        b[IDX_SIGN_LENGTH] = 0;
        b[IDX_EXPONENT] = 0;
        b[IDX_MANTISSA_START] = 0;
        b[IDX_MANTISSA_START + 1] = 0;
        return;
    }
    
    int64_t ma = a[IDX_EXPONENT];
    
    // If exponent is negative (fractional only), result depends on first fraction digit
    if (ma < 0) {
        // |a| < 1, so check if |a| >= 0.5
        if (ma == -1 && a[IDX_MANTISSA_START] >= (RADIX / 2)) {
            // Round to ±1
            b[IDX_PRECISION] = mpnw;
            b[IDX_SIGN_LENGTH] = ia;
            b[IDX_EXPONENT] = 0;
            b[IDX_MANTISSA_START] = 1;
            b[IDX_MANTISSA_START + 1] = 0;
        } else {
            // Round to 0
            b[IDX_PRECISION] = mpnw;
            b[IDX_SIGN_LENGTH] = 0;
            b[IDX_EXPONENT] = 0;
            b[IDX_MANTISSA_START] = 0;
            b[IDX_MANTISSA_START + 1] = 0;
        }
        return;
    }
    
    // Integer part extraction
    int ni = static_cast<int>(ma) + 1;  // Number of integer words
    
    if (ni >= na) {
        // Already an integer (no fractional part)
        mpeq<WORDS>(a, b, mpnw);
        return;
    }
    
    // Copy integer part
    b[IDX_PRECISION] = mpnw;
    b[IDX_EXPONENT] = ma;
    
    for (int i = 0; i < ni; ++i) {
        b[IDX_MANTISSA_START + i] = a[IDX_MANTISSA_START + i];
    }
    
    // Check if we need to round up
    bool round_up = false;
    if (ni < na) {
        // Check first fractional word
        round_up = (a[IDX_MANTISSA_START + ni] >= (RADIX / 2));
    }
    
    // Zero remaining words
    for (int i = ni; i < mpnw + 2; ++i) {
        b[IDX_MANTISSA_START + i] = 0;
    }
    
    b[IDX_SIGN_LENGTH] = combine_sign_length(ia, ni);
    
    if (round_up) {
        // Add 1 to magnitude
        b[IDX_MANTISSA_START + ni - 1] += 1;
        
        // Propagate carry
        for (int i = ni - 1; i >= 0; --i) {
            if (b[IDX_MANTISSA_START + i] >= RADIX) {
                b[IDX_MANTISSA_START + i] -= RADIX;
                if (i > 0) {
                    b[IDX_MANTISSA_START + i - 1] += 1;
                } else {
                    // Overflow - shift right and increment exponent
                    for (int j = ni; j >= 1; --j) {
                        b[IDX_MANTISSA_START + j] = b[IDX_MANTISSA_START + j - 1];
                    }
                    b[IDX_MANTISSA_START] = 1;
                    b[IDX_EXPONENT] += 1;
                    ni = (ni + 1 < mpnw) ? ni + 1 : mpnw;
                    b[IDX_SIGN_LENGTH] = combine_sign_length(ia, ni);
                }
            } else {
                break;
            }
        }
    }
    
    round<WORDS>(b, mpnw);
}

// =============================================================================
// Internal sincos implementation (mpcssnr)
// =============================================================================

/**
 * @brief Compute sine and cosine together
 * 
 * Algorithm (from MPFUN mpcssn):
 * 1. Argument reduction: x = x mod 2π, then to [-π, π]
 * 2. Further reduce: x = x/2^nq (where nq = sqrt(0.5 * precision_bits))
 * 3. Taylor series for sin(x):
 *    sin(s) = s - s³/3! + s⁵/5! - s⁷/7! + ...
 * 4. Double-angle formulas nq times:
 *    First iteration: cos(2x) = 1 - 2*sin²(x)
 *    Subsequent: cos(2x) = 2*cos²(x) - 1
 * 5. sin(t) = ±√(1 - cos²(t))
 * 
 * @tparam WORDS Precision in mantissa words
 * @param a_in Input angle in radians
 * @param cos_out Output: cos(a)
 * @param sin_out Output: sin(a)
 * @param mpnw Working precision
 */
template <int WORDS>
KOKKOS_INLINE_FUNCTION
void mpcssnr(const int64_t* a_in, int64_t* cos_out, int64_t* sin_out, int mpnw) {
    constexpr int ITRMX = 1000000;
    
    int na = extract_length(a_in[IDX_SIGN_LENGTH]);
    if (na > mpnw) na = mpnw;
    
    // Handle zero input: cos(0) = 1, sin(0) = 0
    if (na == 0) {
        // cos = 1
        cos_out[IDX_PRECISION] = mpnw;
        cos_out[IDX_SIGN_LENGTH] = 1;
        cos_out[IDX_EXPONENT] = 0;
        cos_out[IDX_MANTISSA_START] = 1;
        cos_out[IDX_MANTISSA_START + 1] = 0;
        
        // sin = 0
        sin_out[IDX_PRECISION] = mpnw;
        sin_out[IDX_SIGN_LENGTH] = 0;
        sin_out[IDX_EXPONENT] = 0;
        sin_out[IDX_MANTISSA_START] = 0;
        sin_out[IDX_MANTISSA_START + 1] = 0;
        return;
    }
    
    int mpnw1 = mpnw + 1;
    
    // Temporary arrays
    int64_t s0[WORDS + 10];
    int64_t s1[WORDS + 10];
    int64_t s2[WORDS + 10];
    int64_t s3[WORDS + 10];
    int64_t s4[WORDS + 10];
    int64_t s5[WORDS + 10];
    int64_t f1[10];  // Constant 1
    int64_t f2[10];  // Constant 0.5
    int64_t pi_arr[WORDS + 10];
    int64_t two_pi_arr[WORDS + 10];
    
    for (int i = 0; i < WORDS + 10; ++i) {
        s0[i] = s1[i] = s2[i] = s3[i] = s4[i] = s5[i] = pi_arr[i] = two_pi_arr[i] = 0;
    }
    for (int i = 0; i < 10; ++i) {
        f1[i] = f2[i] = 0;
    }
    
    s0[IDX_ALLOCATED] = mpnw + 7;
    s1[IDX_ALLOCATED] = mpnw + 7;
    s2[IDX_ALLOCATED] = mpnw + 7;
    s3[IDX_ALLOCATED] = mpnw + 7;
    s4[IDX_ALLOCATED] = mpnw + 7;
    s5[IDX_ALLOCATED] = mpnw + 7;
    pi_arr[IDX_ALLOCATED] = mpnw + 7;
    two_pi_arr[IDX_ALLOCATED] = mpnw + 7;
    
    // Set f1 = 1
    f1[IDX_ALLOCATED] = 9;
    f1[IDX_PRECISION] = mpnw;
    f1[IDX_SIGN_LENGTH] = 1;
    f1[IDX_EXPONENT] = 0;
    f1[IDX_MANTISSA_START] = 1;
    f1[IDX_MANTISSA_START + 1] = 0;
    
    // Set f2 = 0.5
    f2[IDX_ALLOCATED] = 9;
    f2[IDX_PRECISION] = mpnw;
    f2[IDX_SIGN_LENGTH] = 1;
    f2[IDX_EXPONENT] = -1;
    f2[IDX_MANTISSA_START] = RADIX / 2;  // 0.5 * RADIX in first word with exp=-1
    f2[IDX_MANTISSA_START + 1] = 0;
    
    // Get pi constant - for now, use high-precision approximation
    // In production, this would come from precomputed constants
    // Pi ≈ 3.14159265358979323846...
    mpdmc<WORDS>(3.141592653589793238, 0, pi_arr, mpnw1);
    
    // Compute 2*pi
    muld<WORDS>(pi_arr, 2.0, two_pi_arr, mpnw1);
    
    // Step 1: Reduce argument to [-π, π]
    // s1 = a / (2π), s2 = nint(s1), s3 = a - 2π * s2
    div<WORDS>(a_in, two_pi_arr, s1, mpnw1);
    mpnint<WORDS>(s1, s2, mpnw1);
    mul<WORDS>(two_pi_arr, s2, s4, mpnw1);
    sub<WORDS>(a_in, s4, s3, mpnw1);
    
    // Check if reduced argument is zero
    if (s3[IDX_SIGN_LENGTH] == 0) {
        // cos = 1, sin = 0
        cos_out[IDX_PRECISION] = mpnw;
        cos_out[IDX_SIGN_LENGTH] = 1;
        cos_out[IDX_EXPONENT] = 0;
        cos_out[IDX_MANTISSA_START] = 1;
        cos_out[IDX_MANTISSA_START + 1] = 0;
        
        sin_out[IDX_PRECISION] = mpnw;
        sin_out[IDX_SIGN_LENGTH] = 0;
        sin_out[IDX_EXPONENT] = 0;
        sin_out[IDX_MANTISSA_START] = 0;
        sin_out[IDX_MANTISSA_START + 1] = 0;
        return;
    }
    
    // Step 2: Determine nq for further reduction
    // If argument already very small (exponent < -1), skip division
    int nq;
    if (s3[IDX_EXPONENT] >= -1) {
        // nq = int(sqrt(0.5 * mpnw1 * BITS_PER_WORD))
        double t = 0.5 * mpnw1 * BITS_PER_WORD;
        nq = static_cast<int>(Kokkos::sqrt(t));
        if (nq < 1) nq = 1;
    } else {
        nq = 0;
    }
    
    // Divide by 2^nq
    if (nq > 0) {
        double divisor = Kokkos::pow(2.0, static_cast<double>(nq));
        divd<WORDS>(s3, divisor, s0, mpnw1);
    } else {
        mpeq<WORDS>(s3, s0, mpnw1);
    }
    mpeq<WORDS>(s0, s1, mpnw1);  // s1 = reduced argument s
    
    // Save sign of original reduced sin for later
    int is = extract_sign(s0[IDX_SIGN_LENGTH]);
    
    // Step 3: Taylor series for sin(s)
    // sin(s) = s - s³/3! + s⁵/5! - ...
    mul<WORDS>(s0, s0, s2, mpnw1);  // s2 = s²
    
    int mpnw2 = mpnw1;
    
    for (int i1 = 1; i1 <= ITRMX; ++i1) {
        double t2 = -static_cast<double>(2 * i1) * static_cast<double>(2 * i1 + 1);
        mul<WORDS>(s2, s1, s3, mpnw2);      // s3 = s² * s_{prev}
        divd<WORDS>(s3, t2, s1, mpnw2);     // s1 = s3 / (2k)(2k+1)
        add<WORDS>(s1, s0, s3, mpnw1);      // s3 = s0 + s1
        mpeq<WORDS>(s3, s0, mpnw1);         // s0 = s3
        
        // Check convergence
        if (s1[IDX_SIGN_LENGTH] == 0) break;
        int64_t exp_s1 = s1[IDX_EXPONENT];
        int64_t exp_s0 = s0[IDX_EXPONENT];
        if (exp_s1 < exp_s0 - mpnw1) break;
        
        // Adjust working precision
        int new_mpnw2 = mpnw1 + static_cast<int>(exp_s1 - exp_s0) + 1;
        if (new_mpnw2 < 4) new_mpnw2 = 4;
        if (new_mpnw2 > mpnw1) new_mpnw2 = mpnw1;
        mpnw2 = new_mpnw2;
    }
    
    // s0 now contains sin(s) where s = reduced_arg / 2^nq
    
    // Step 4: Apply double-angle formulas nq times
    if (nq > 0) {
        // First iteration: cos(2x) = 1 - 2*sin²(x) = 2*(0.5 - sin²(x))
        mul<WORDS>(s0, s0, s4, mpnw1);       // s4 = sin²(x)
        sub<WORDS>(f2, s4, s5, mpnw1);       // s5 = 0.5 - sin²(x)
        muld<WORDS>(s5, 2.0, s0, mpnw1);     // s0 = cos(2x)
        
        // Subsequent iterations: cos(2x) = 2*cos²(x) - 1 = 2*(cos²(x) - 0.5)
        for (int j = 2; j <= nq; ++j) {
            mul<WORDS>(s0, s0, s4, mpnw1);   // s4 = cos²(x)
            sub<WORDS>(s4, f2, s5, mpnw1);   // s5 = cos²(x) - 0.5
            muld<WORDS>(s5, 2.0, s0, mpnw1); // s0 = cos(2x)
        }
        
        // Now s0 = cos(reduced_arg)
        // Compute sin = ±√(1 - cos²)
        mul<WORDS>(s0, s0, s4, mpnw1);       // s4 = cos²
        sub<WORDS>(f1, s4, s5, mpnw1);       // s5 = 1 - cos²
        sqrt<WORDS>(s5, s1, mpnw1);          // s1 = |sin|
        
        // Correct sign of sin
        if (is < 0) {
            s1[IDX_SIGN_LENGTH] = -extract_length(s1[IDX_SIGN_LENGTH]);
        }
    } else {
        // nq = 0: s0 = sin(reduced_arg)
        mpeq<WORDS>(s0, s1, mpnw1);          // s1 = sin
        
        // Compute cos = √(1 - sin²)
        mul<WORDS>(s0, s0, s4, mpnw1);       // s4 = sin²
        sub<WORDS>(f1, s4, s5, mpnw1);       // s5 = 1 - sin²
        sqrt<WORDS>(s5, s0, mpnw1);          // s0 = cos
    }
    
    // Round and copy results
    round<WORDS>(s0, mpnw);
    round<WORDS>(s1, mpnw);
    mpeq<WORDS>(s0, cos_out, mpnw);
    mpeq<WORDS>(s1, sin_out, mpnw);
}

// =============================================================================
// Internal atan2 implementation (mpang)
// =============================================================================

/**
 * @brief Compute atan2(y, x) - the angle for point (x, y)
 * 
 * Algorithm (from MPFUN mpang):
 * 1. Handle special cases (x=0, y=0)
 * 2. Normalize: x' = x/r, y' = y/r where r = √(x² + y²)
 * 3. Initial approximation using double precision
 * 4. Newton iteration:
 *    If |x| <= |y|: z_{k+1} = z_k - [x - cos(z_k)] / sin(z_k)
 *    If |y| < |x|:  z_{k+1} = z_k + [y - sin(z_k)] / cos(z_k)
 * 
 * @tparam WORDS Precision in mantissa words
 * @param x_in X coordinate
 * @param y_in Y coordinate
 * @param a_out Output angle in radians (-π, π]
 * @param mpnw Working precision
 */
template <int WORDS>
KOKKOS_INLINE_FUNCTION
void mpang(const int64_t* x_in, const int64_t* y_in, int64_t* a_out, int mpnw) {
    constexpr int NIT = 3;  // Extra iterations at full precision
    
    int ix = extract_sign(x_in[IDX_SIGN_LENGTH]);
    int iy = extract_sign(y_in[IDX_SIGN_LENGTH]);
    int nx = extract_length(x_in[IDX_SIGN_LENGTH]);
    int ny = extract_length(y_in[IDX_SIGN_LENGTH]);
    
    if (nx > mpnw) nx = mpnw;
    if (ny > mpnw) ny = mpnw;
    
    int mpnw1 = mpnw + 1;
    
    // Both x and y are zero - undefined
    if (nx == 0 && ny == 0) {
        a_out[IDX_ALLOCATED] = -1;  // NaN
        a_out[IDX_PRECISION] = mpnw;
        a_out[IDX_SIGN_LENGTH] = 0;
        return;
    }
    
    // Temporary arrays
    int64_t s0[WORDS + 10];
    int64_t s1[WORDS + 10];
    int64_t s2[WORDS + 10];
    int64_t s3[WORDS + 10];
    int64_t s4[WORDS + 10];
    int64_t s5[WORDS + 10];
    int64_t pi_arr[WORDS + 10];
    
    for (int i = 0; i < WORDS + 10; ++i) {
        s0[i] = s1[i] = s2[i] = s3[i] = s4[i] = s5[i] = pi_arr[i] = 0;
    }
    
    s0[IDX_ALLOCATED] = mpnw + 7;
    s1[IDX_ALLOCATED] = mpnw + 7;
    s2[IDX_ALLOCATED] = mpnw + 7;
    s3[IDX_ALLOCATED] = mpnw + 7;
    s4[IDX_ALLOCATED] = mpnw + 7;
    s5[IDX_ALLOCATED] = mpnw + 7;
    pi_arr[IDX_ALLOCATED] = mpnw + 7;
    
    // Get pi constant
    mpdmc<WORDS>(3.141592653589793238, 0, pi_arr, mpnw1);
    
    // Special case: x = 0
    if (nx == 0) {
        if (iy > 0) {
            // +π/2
            muld<WORDS>(pi_arr, 0.5, a_out, mpnw);
        } else {
            // -π/2
            muld<WORDS>(pi_arr, -0.5, a_out, mpnw);
        }
        return;
    }
    
    // Special case: y = 0
    if (ny == 0) {
        if (ix > 0) {
            // 0
            a_out[IDX_PRECISION] = mpnw;
            a_out[IDX_SIGN_LENGTH] = 0;
            a_out[IDX_EXPONENT] = 0;
            a_out[IDX_MANTISSA_START] = 0;
            a_out[IDX_MANTISSA_START + 1] = 0;
        } else {
            // π
            mpeq<WORDS>(pi_arr, a_out, mpnw);
        }
        return;
    }
    
    // Determine MQ = ceil(log2(mpnw))
    double t1 = static_cast<double>(mpnw1);
    int mq = static_cast<int>(LOG2_E * Kokkos::log(t1) + 2.0 - COMPARE_FUZZ);
    
    // Normalize x and y: compute r = √(x² + y²)
    mul<WORDS>(x_in, x_in, s0, mpnw1);  // s0 = x²
    mul<WORDS>(y_in, y_in, s1, mpnw1);  // s1 = y²
    add<WORDS>(s0, s1, s2, mpnw1);      // s2 = x² + y²
    sqrt<WORDS>(s2, s3, mpnw1);         // s3 = r = √(x² + y²)
    div<WORDS>(x_in, s3, s1, mpnw1);    // s1 = x/r = cos(angle)
    div<WORDS>(y_in, s3, s2, mpnw1);    // s2 = y/r = sin(angle)
    
    // Initial approximation using double precision
    double d1, d2;
    int n1, n2;
    mpmdc<WORDS>(s1, d1, n1, mpnw1);
    mpmdc<WORDS>(s2, d2, n2, mpnw1);
    
    // Clamp exponents
    if (n1 < -BITS_PER_WORD) n1 = -BITS_PER_WORD;
    if (n2 < -BITS_PER_WORD) n2 = -BITS_PER_WORD;
    
    d1 = d1 * Kokkos::pow(2.0, static_cast<double>(n1));
    d2 = d2 * Kokkos::pow(2.0, static_cast<double>(n2));
    
    double t3 = Kokkos::atan2(d2, d1);
    mpdmc<WORDS>(t3, 0, s5, mpnw1);  // s5 = initial angle approximation
    
    // Choose Newton iteration based on which has larger denominator
    int kk;
    if (Kokkos::fabs(d1) <= Kokkos::fabs(d2)) {
        kk = 1;  // Use cos iteration: z -= (x - cos(z)) / sin(z)
        mpeq<WORDS>(s1, s0, mpnw1);  // s0 = target (x/r = cos)
    } else {
        kk = 2;  // Use sin iteration: z += (y - sin(z)) / cos(z)
        mpeq<WORDS>(s2, s0, mpnw1);  // s0 = target (y/r = sin)
    }
    
    // Newton-Raphson iteration with dynamically increasing precision
    int mpnw_iter = 4;
    int iq = 0;
    
    for (int k = 1; k <= mq; ++k) {
        mpnw_iter = (2 * mpnw_iter - 2 < mpnw) ? 2 * mpnw_iter - 2 : mpnw;
        if (mpnw_iter < 4) mpnw_iter = 4;
        
        do {
            // Compute sin and cos of current angle estimate
            mpcssnr<WORDS>(s5, s1, s2, mpnw_iter);  // s1 = cos(z), s2 = sin(z)
            
            if (kk == 1) {
                // z -= (x - cos(z)) / sin(z)
                sub<WORDS>(s0, s1, s3, mpnw_iter);  // s3 = x - cos(z)
                div<WORDS>(s3, s2, s4, mpnw_iter);  // s4 = (x - cos(z)) / sin(z)
                sub<WORDS>(s5, s4, s1, mpnw_iter);  // s1 = z - correction
            } else {
                // z += (y - sin(z)) / cos(z)
                sub<WORDS>(s0, s2, s3, mpnw_iter);  // s3 = y - sin(z)
                div<WORDS>(s3, s1, s4, mpnw_iter);  // s4 = (y - sin(z)) / cos(z)
                add<WORDS>(s5, s4, s1, mpnw_iter);  // s1 = z + correction
            }
            mpeq<WORDS>(s1, s5, mpnw_iter);
            
            if (k == mq - NIT && iq == 0) {
                iq = 1;
            } else {
                break;
            }
        } while (true);
    }
    
    // Round and copy result
    round<WORDS>(s5, mpnw);
    mpeq<WORDS>(s5, a_out, mpnw);
}

} // namespace detail

// =============================================================================
// Public API: sincos(x, &sin_out, &cos_out)
// =============================================================================

/**
 * @brief Simultaneous sine and cosine computation
 * 
 * More efficient than separate sin() and cos() calls as it shares
 * the argument reduction work.
 * 
 * @tparam W Number of mantissa words (precision)
 * @param x Angle in radians
 * @param sin_out Output: sin(x)
 * @param cos_out Output: cos(x)
 */
template <int W>
KOKKOS_INLINE_FUNCTION 
void sincos(const MPFloat<W>& x, MPFloat<W>& sin_out, MPFloat<W>& cos_out) {
    detail::mpcssnr<W>(x.data(), cos_out.data(), sin_out.data(), W);
}

// =============================================================================
// Public API: sin(x)
// =============================================================================

/**
 * @brief Sine function
 * 
 * @tparam W Number of mantissa words (precision)
 * @param x Angle in radians
 * @return sin(x)
 */
template <int W>
KOKKOS_INLINE_FUNCTION 
MPFloat<W> sin(const MPFloat<W>& x) {
    MPFloat<W> sin_out, cos_out;
    sincos(x, sin_out, cos_out);
    return sin_out;
}

// =============================================================================
// Public API: cos(x)
// =============================================================================

/**
 * @brief Cosine function
 * 
 * @tparam W Number of mantissa words (precision)
 * @param x Angle in radians
 * @return cos(x)
 */
template <int W>
KOKKOS_INLINE_FUNCTION 
MPFloat<W> cos(const MPFloat<W>& x) {
    MPFloat<W> sin_out, cos_out;
    sincos(x, sin_out, cos_out);
    return cos_out;
}

// =============================================================================
// Public API: tan(x)
// =============================================================================

/**
 * @brief Tangent function
 * 
 * @tparam W Number of mantissa words (precision)
 * @param x Angle in radians
 * @return tan(x) = sin(x) / cos(x)
 * @note Returns NaN when cos(x) ≈ 0
 */
template <int W>
KOKKOS_INLINE_FUNCTION 
MPFloat<W> tan(const MPFloat<W>& x) {
    MPFloat<W> sin_out, cos_out;
    sincos(x, sin_out, cos_out);
    return sin_out / cos_out;
}

// =============================================================================
// Public API: atan(x)
// =============================================================================

/**
 * @brief Arctangent function
 * 
 * @tparam W Number of mantissa words (precision)
 * @param x Input value (any real)
 * @return arctan(x) in [-π/2, π/2]
 */
template <int W>
KOKKOS_INLINE_FUNCTION 
MPFloat<W> atan(const MPFloat<W>& x) {
    // atan(x) = atan2(x, 1)
    MPFloat<W> one(1.0);
    MPFloat<W> result;
    detail::mpang<W>(one.data(), x.data(), result.data(), W);
    return result;
}

// =============================================================================
// Public API: atan2(y, x)
// =============================================================================

/**
 * @brief Two-argument arctangent
 * 
 * @tparam W Number of mantissa words (precision)
 * @param y The y-coordinate (numerator)
 * @param x The x-coordinate (denominator)
 * @return arctan(y/x) with correct quadrant, result in [-π, π]
 */
template <int W>
KOKKOS_INLINE_FUNCTION 
MPFloat<W> atan2(const MPFloat<W>& y, const MPFloat<W>& x) {
    MPFloat<W> result;
    detail::mpang<W>(x.data(), y.data(), result.data(), W);
    return result;
}

// =============================================================================
// Public API: asin(x)
// =============================================================================

/**
 * @brief Arcsine function
 * 
 * @tparam W Number of mantissa words (precision)
 * @param x Input value in [-1, 1]
 * @return arcsin(x) in [-π/2, π/2]
 * @pre |x| <= 1 (returns NaN otherwise)
 * @note Uses identity: asin(x) = atan(x / √(1-x²))
 */
template <int W>
KOKKOS_INLINE_FUNCTION 
MPFloat<W> asin(const MPFloat<W>& x) {
    // asin(x) = atan2(x, √(1 - x²))
    MPFloat<W> one(1.0);
    MPFloat<W> x_sq = x * x;
    MPFloat<W> one_minus_x_sq = one - x_sq;
    
    // Check for |x| > 1
    if (one_minus_x_sq.is_negative()) {
        MPFloat<W> result;
        result.set_nan();
        return result;
    }
    
    MPFloat<W> sqrt_term = sqrt(one_minus_x_sq);
    return atan2(x, sqrt_term);
}

// =============================================================================
// Public API: acos(x)
// =============================================================================

/**
 * @brief Arccosine function
 * 
 * @tparam W Number of mantissa words (precision)
 * @param x Input value in [-1, 1]
 * @return arccos(x) in [0, π]
 * @pre |x| <= 1 (returns NaN otherwise)
 * @note Uses identity: acos(x) = atan2(√(1-x²), x)
 */
template <int W>
KOKKOS_INLINE_FUNCTION 
MPFloat<W> acos(const MPFloat<W>& x) {
    // acos(x) = atan2(√(1 - x²), x)
    MPFloat<W> one(1.0);
    MPFloat<W> x_sq = x * x;
    MPFloat<W> one_minus_x_sq = one - x_sq;
    
    // Check for |x| > 1
    if (one_minus_x_sq.is_negative()) {
        MPFloat<W> result;
        result.set_nan();
        return result;
    }
    
    MPFloat<W> sqrt_term = sqrt(one_minus_x_sq);
    return atan2(sqrt_term, x);
}

} // namespace mpfun

#endif // MPFUN_TRANSCENDENTAL_TRIG_HPP
