/**
 * @file mp_float.hpp
 * @brief Core MPFloat<WORDS> arbitrary precision floating-point type
 * 
 * This file defines the MPFloat class template, which provides fixed-size
 * arbitrary precision floating-point arithmetic compatible with Kokkos
 * parallel execution on CPUs and GPUs.
 * 
 * Template parameter WORDS determines the precision:
 * - Each word provides ~18.06 decimal digits
 * - MPFloat<6> provides ~100 decimal digits
 * - MPFloat<56> provides ~1000 decimal digits
 * 
 * Copyright (c) 2024 Argonne National Laboratory
 * Ported from MPFUN2020-Fort by David H. Bailey
 */

#ifndef MPFUN_MP_FLOAT_HPP
#define MPFUN_MP_FLOAT_HPP

#include "core/representation.hpp"
#include <cstdint>
#include <cmath>

// Kokkos mock functions are defined in core/representation.hpp
// No additional mock functions needed here

namespace mpfun {

// Forward declarations
template <int WORDS> class MPFloat;

template <int W>
KOKKOS_INLINE_FUNCTION MPFloat<W> abs(const MPFloat<W>& x);

template <int W>
KOKKOS_INLINE_FUNCTION int compare(const MPFloat<W>& a, const MPFloat<W>& b);

/**
 * @brief Fixed-precision arbitrary-precision floating-point number
 * 
 * MPFloat<WORDS> stores a number as:
 *   value = sign * mantissa * RADIX^exponent
 * 
 * where RADIX = 2^60 and mantissa is represented as an array of WORDS
 * 60-bit integers.
 * 
 * @tparam WORDS Number of mantissa words (determines precision)
 */
template <int WORDS>
class MPFloat {
public:
    static_assert(WORDS >= 1, "MPFloat requires at least 1 word");
    static_assert(WORDS <= 200000, "MPFloat word count exceeds maximum");
    
    // Type aliases
    using word_type = int64_t;
    
    //=========================================================================
    // Compile-time constants
    //=========================================================================
    
    /// Number of mantissa words
    static constexpr int capacity() { return WORDS; }
    
    /// Approximate number of decimal digits of precision
    static constexpr int approx_digits() {
        return static_cast<int>(WORDS * detail::DIGITS_PER_WORD);
    }
    
    /// Total size of the data array
    /// MPFUN operations require mpnw+6 or mpnw+7 elements, so we allocate
    /// METADATA_WORDS + WORDS + 3 = WORDS + 7 to provide adequate padding.
    static constexpr int data_size() { return WORDS + detail::METADATA_WORDS + 3; }
    
    //=========================================================================
    // Constructors
    //=========================================================================
    
    /**
     * @brief Default constructor - initializes to zero
     */
    KOKKOS_INLINE_FUNCTION
    MPFloat() {
        data_[detail::IDX_ALLOCATED] = WORDS;
        data_[detail::IDX_PRECISION] = WORDS;
        data_[detail::IDX_SIGN_LENGTH] = 0;
        data_[detail::IDX_EXPONENT] = 0;
        data_[detail::IDX_MANTISSA_START] = 0;
        data_[detail::IDX_MANTISSA_START + 1] = 0;
    }
    
    /**
     * @brief Construct from double precision value
     * @param d Double precision value to convert
     */
    KOKKOS_INLINE_FUNCTION
    explicit MPFloat(double d) {
        from_double(d, 0);
    }
    
    /**
     * @brief Construct from double with power-of-2 scaling
     * @param d Double precision value
     * @param n Power of 2 to multiply by (result = d * 2^n)
     */
    KOKKOS_INLINE_FUNCTION
    MPFloat(double d, int n) {
        from_double(d, n);
    }
    
    /**
     * @brief Construct from 64-bit integer
     * @param i Integer value to convert
     */
    KOKKOS_INLINE_FUNCTION
    explicit MPFloat(int64_t i) {
        from_double(static_cast<double>(i), 0);
    }
    
    /**
     * @brief Construct from 32-bit integer
     * @param i Integer value to convert
     */
    KOKKOS_INLINE_FUNCTION
    explicit MPFloat(int i) {
        from_double(static_cast<double>(i), 0);
    }
    
    /**
     * @brief Copy constructor
     */
    KOKKOS_INLINE_FUNCTION
    MPFloat(const MPFloat& other) = default;
    
    /**
     * @brief Copy assignment
     */
    KOKKOS_INLINE_FUNCTION
    MPFloat& operator=(const MPFloat& other) = default;
    
    //=========================================================================
    // Assignment operators
    //=========================================================================
    
    /**
     * @brief Assign from double
     */
    KOKKOS_INLINE_FUNCTION
    MPFloat& operator=(double d) {
        from_double(d, 0);
        return *this;
    }
    
    //=========================================================================
    // Accessor methods
    //=========================================================================
    
    /**
     * @brief Get the sign of the number
     * @return -1, 0, or +1
     */
    KOKKOS_INLINE_FUNCTION
    int sign() const {
        return detail::extract_sign(data_[detail::IDX_SIGN_LENGTH]);
    }
    
    /**
     * @brief Get the exponent (power of RADIX)
     */
    KOKKOS_INLINE_FUNCTION
    int64_t exponent() const {
        return data_[detail::IDX_EXPONENT];
    }
    
    /**
     * @brief Get the number of mantissa words currently used
     */
    KOKKOS_INLINE_FUNCTION
    int mantissa_length() const {
        return detail::extract_length(data_[detail::IDX_SIGN_LENGTH]);
    }
    
    /**
     * @brief Get the working precision in words
     */
    KOKKOS_INLINE_FUNCTION
    int precision_words() const {
        return static_cast<int>(data_[detail::IDX_PRECISION]);
    }
    
    /**
     * @brief Check if the value is zero
     */
    KOKKOS_INLINE_FUNCTION
    bool is_zero() const {
        return data_[detail::IDX_SIGN_LENGTH] == 0;
    }
    
    /**
     * @brief Check if the value is negative
     */
    KOKKOS_INLINE_FUNCTION
    bool is_negative() const {
        return data_[detail::IDX_SIGN_LENGTH] < 0;
    }
    
    /**
     * @brief Check if the value is positive
     */
    KOKKOS_INLINE_FUNCTION
    bool is_positive() const {
        return data_[detail::IDX_SIGN_LENGTH] > 0;
    }
    
    /**
     * @brief Check if this represents a NaN (invalid state)
     */
    KOKKOS_INLINE_FUNCTION
    bool is_nan() const {
        return data_[detail::IDX_ALLOCATED] < 0;
    }
    
    /**
     * @brief Set this value to NaN (invalid state)
     */
    KOKKOS_INLINE_FUNCTION
    void set_nan() {
        data_[detail::IDX_ALLOCATED] = -1;
    }
    
    //=========================================================================
    // Conversion methods
    //=========================================================================
    
    /**
     * @brief Convert to double precision approximation
     * @return Double precision approximation of this value
     */
    KOKKOS_INLINE_FUNCTION
    double to_double() const {
        if (is_zero()) return 0.0;
        
        int na = mantissa_length();
        double aa = static_cast<double>(data_[detail::IDX_MANTISSA_START]);
        if (na >= 2) {
            aa += static_cast<double>(data_[detail::IDX_MANTISSA_START + 1]) / 
                  static_cast<double>(detail::RADIX);
        }
        
        int64_t exp = data_[detail::IDX_EXPONENT];
        int n = static_cast<int>(detail::BITS_PER_WORD * exp);
        
        // Apply sign
        aa = sign() > 0 ? aa : -aa;
        
        // Reduce to within [1, 2)
        int shift = static_cast<int>(Kokkos::log(Kokkos::fabs(aa)) / Kokkos::log(2.0) + 
                                      detail::COMPARE_FUZZ);
        aa = aa / Kokkos::pow(2.0, static_cast<double>(shift));
        n += shift;
        
        if (Kokkos::fabs(aa) < 1.0) {
            aa *= 2.0;
            n -= 1;
        } else if (Kokkos::fabs(aa) >= 2.0) {
            aa *= 0.5;
            n += 1;
        }
        
        return aa * Kokkos::pow(2.0, static_cast<double>(n));
    }
    
    //=========================================================================
    // Arithmetic operators (stub declarations - implemented in core headers)
    //=========================================================================
    
    KOKKOS_INLINE_FUNCTION
    MPFloat operator+(const MPFloat& rhs) const;
    
    KOKKOS_INLINE_FUNCTION
    MPFloat operator-(const MPFloat& rhs) const;
    
    KOKKOS_INLINE_FUNCTION
    MPFloat operator*(const MPFloat& rhs) const;
    
    KOKKOS_INLINE_FUNCTION
    MPFloat operator/(const MPFloat& rhs) const;
    
    /**
     * @brief Unary negation
     */
    KOKKOS_INLINE_FUNCTION
    MPFloat operator-() const {
        MPFloat result(*this);
        int na = mantissa_length();
        int s = sign();
        result.data_[detail::IDX_SIGN_LENGTH] = detail::combine_sign_length(-s, na);
        return result;
    }
    
    //=========================================================================
    // Compound assignment operators
    //=========================================================================
    
    KOKKOS_INLINE_FUNCTION
    MPFloat& operator+=(const MPFloat& rhs) {
        *this = *this + rhs;
        return *this;
    }
    
    KOKKOS_INLINE_FUNCTION
    MPFloat& operator-=(const MPFloat& rhs) {
        *this = *this - rhs;
        return *this;
    }
    
    KOKKOS_INLINE_FUNCTION
    MPFloat& operator*=(const MPFloat& rhs) {
        *this = *this * rhs;
        return *this;
    }
    
    KOKKOS_INLINE_FUNCTION
    MPFloat& operator/=(const MPFloat& rhs) {
        *this = *this / rhs;
        return *this;
    }
    
    //=========================================================================
    // Comparison operators
    //=========================================================================
    
    KOKKOS_INLINE_FUNCTION
    bool operator==(const MPFloat& rhs) const {
        return compare(*this, rhs) == 0;
    }
    
    KOKKOS_INLINE_FUNCTION
    bool operator!=(const MPFloat& rhs) const {
        return compare(*this, rhs) != 0;
    }
    
    KOKKOS_INLINE_FUNCTION
    bool operator<(const MPFloat& rhs) const {
        return compare(*this, rhs) < 0;
    }
    
    KOKKOS_INLINE_FUNCTION
    bool operator<=(const MPFloat& rhs) const {
        return compare(*this, rhs) <= 0;
    }
    
    KOKKOS_INLINE_FUNCTION
    bool operator>(const MPFloat& rhs) const {
        return compare(*this, rhs) > 0;
    }
    
    KOKKOS_INLINE_FUNCTION
    bool operator>=(const MPFloat& rhs) const {
        return compare(*this, rhs) >= 0;
    }
    
    //=========================================================================
    // Raw data access (use with caution)
    //=========================================================================
    
    KOKKOS_INLINE_FUNCTION
    word_type* data() { return data_; }
    
    KOKKOS_INLINE_FUNCTION
    const word_type* data() const { return data_; }
    
    //=========================================================================
    // Friend declarations
    //=========================================================================
    
    template <int W>
    friend KOKKOS_INLINE_FUNCTION MPFloat<W> abs(const MPFloat<W>& x);
    
    template <int W>
    friend KOKKOS_INLINE_FUNCTION int compare(const MPFloat<W>& a, const MPFloat<W>& b);
    
private:
    /// The data array: metadata + mantissa + guard words
    /// MPFUN operations require mpnw+6 or mpnw+7 elements, so we need extra space
    word_type data_[WORDS + detail::METADATA_WORDS + 3];
    
    //=========================================================================
    // Internal helper methods
    //=========================================================================
    
    /**
     * @brief Convert double to internal representation
     * 
     * Implements mpdmc from mpfunb.f90: converts a * 2^n to MPFloat.
     * 
     * @param a Double precision value
     * @param n Power of 2 scaling
     */
    KOKKOS_INLINE_FUNCTION
    void from_double(double a, int n) {
        data_[detail::IDX_ALLOCATED] = WORDS;
        data_[detail::IDX_PRECISION] = WORDS;
        
        // Check for zero
        if (a == 0.0) {
            data_[detail::IDX_SIGN_LENGTH] = 0;
            data_[detail::IDX_EXPONENT] = 0;
            data_[detail::IDX_MANTISSA_START] = 0;
            data_[detail::IDX_MANTISSA_START + 1] = 0;
            return;
        }
        
        int n1 = n / detail::BITS_PER_WORD;
        int n2 = n - detail::BITS_PER_WORD * n1;
        double aa = Kokkos::fabs(a) * Kokkos::pow(2.0, static_cast<double>(n2));
        
        // Reduce aa to within [1, RADIX)
        const double bdx = static_cast<double>(detail::RADIX);
        
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
        
        // Store the result
        data_[detail::IDX_EXPONENT] = n1;
        data_[detail::IDX_MANTISSA_START] = static_cast<int64_t>(aa);
        aa = bdx * (aa - static_cast<double>(data_[detail::IDX_MANTISSA_START]));
        data_[detail::IDX_MANTISSA_START + 1] = static_cast<int64_t>(aa);
        data_[detail::IDX_MANTISSA_START + 2] = 0;
        data_[detail::IDX_MANTISSA_START + 3] = 0;
        data_[detail::IDX_MANTISSA_START + 4] = 0;
        
        // Determine actual length
        int len = 1;
        if (data_[detail::IDX_MANTISSA_START + 1] != 0) len = 2;
        
        data_[detail::IDX_SIGN_LENGTH] = (a >= 0.0) ? len : -len;
    }
    
    /**
     * @brief Set this value to zero
     */
    KOKKOS_INLINE_FUNCTION
    void set_zero() {
        data_[detail::IDX_SIGN_LENGTH] = 0;
        data_[detail::IDX_EXPONENT] = 0;
        data_[detail::IDX_MANTISSA_START] = 0;
        data_[detail::IDX_MANTISSA_START + 1] = 0;
    }
    
    /**
     * @brief Normalize the number after arithmetic operations
     * 
     * This handles carry propagation, leading zero removal, and rounding.
     * Implements mpnorm from mpfunb.f90.
     */
    KOKKOS_INLINE_FUNCTION
    void normalize(int mpnw);
};

//=============================================================================
// Free Functions
//=============================================================================

/**
 * @brief Absolute value
 */
template <int W>
KOKKOS_INLINE_FUNCTION
MPFloat<W> abs(const MPFloat<W>& x) {
    MPFloat<W> result(x);
    int len = x.mantissa_length();
    result.data_[detail::IDX_SIGN_LENGTH] = len;  // Always positive
    return result;
}

/**
 * @brief Maximum of two values
 */
template <int W>
KOKKOS_INLINE_FUNCTION
MPFloat<W> max(const MPFloat<W>& a, const MPFloat<W>& b) {
    return (a > b) ? a : b;
}

/**
 * @brief Minimum of two values
 */
template <int W>
KOKKOS_INLINE_FUNCTION
MPFloat<W> min(const MPFloat<W>& a, const MPFloat<W>& b) {
    return (a < b) ? a : b;
}

//=============================================================================
// Type Aliases for Common Precision Levels
//=============================================================================

/// ~50 decimal digits
using MPFloat50   = MPFloat<3>;

/// ~100 decimal digits  
using MPFloat100  = MPFloat<6>;

/// ~200 decimal digits
using MPFloat200  = MPFloat<12>;

/// ~500 decimal digits
using MPFloat500  = MPFloat<28>;

/// ~1000 decimal digits
using MPFloat1000 = MPFloat<56>;

// Short aliases
using mp50   = MPFloat50;
using mp100  = MPFloat100;
using mp200  = MPFloat200;
using mp500  = MPFloat500;
using mp1000 = MPFloat1000;

} // namespace mpfun

// Include implementation headers
#include "core/add.hpp"
#include "core/mul.hpp"
#include "core/div.hpp"
#include "core/sqrt.hpp"
#include "core/compare.hpp"

//=============================================================================
// Remaining operator implementations (after all includes)
//=============================================================================

namespace mpfun {

template <int WORDS>
KOKKOS_INLINE_FUNCTION
MPFloat<WORDS> MPFloat<WORDS>::operator+(const MPFloat<WORDS>& rhs) const {
    MPFloat<WORDS> result;
    detail::add<WORDS>(data_, rhs.data_, result.data_, WORDS);
    return result;
}

template <int WORDS>
KOKKOS_INLINE_FUNCTION
MPFloat<WORDS> MPFloat<WORDS>::operator-(const MPFloat<WORDS>& rhs) const {
    MPFloat<WORDS> result;
    detail::sub<WORDS>(data_, rhs.data_, result.data_, WORDS);
    return result;
}

template <int WORDS>
KOKKOS_INLINE_FUNCTION
MPFloat<WORDS> MPFloat<WORDS>::operator*(const MPFloat<WORDS>& rhs) const {
    MPFloat<WORDS> result;
    detail::mul<WORDS>(data_, rhs.data_, result.data_, WORDS);
    return result;
}

template <int WORDS>
KOKKOS_INLINE_FUNCTION
MPFloat<WORDS> MPFloat<WORDS>::operator/(const MPFloat<WORDS>& rhs) const {
    MPFloat<WORDS> result;
    detail::div<WORDS>(data_, rhs.data_, result.data_, WORDS);
    return result;
}

} // namespace mpfun

#endif // MPFUN_MP_FLOAT_HPP
