/**
 * @file representation.hpp
 * @brief Memory layout and constants for MPFloat arbitrary precision numbers
 * 
 * Memory Layout (following MPFUN20-Fort):
 * ----------------------------------------
 * data_[0] = allocated_words
 * data_[1] = working_precision
 * data_[2] = sign * length
 * data_[3] = exponent
 * data_[4+] = mantissa words (most significant first)
 * 
 * Copyright (c) 2024 Argonne National Laboratory
 * Ported from MPFUN2020-Fort by David H. Bailey
 */

#ifndef MPFUN_CORE_REPRESENTATION_HPP
#define MPFUN_CORE_REPRESENTATION_HPP

#include <cstdint>
#include <cmath>

// KOKKOS_INLINE_FUNCTION is defined by Kokkos; for standalone testing
// define it if not already present
#ifndef KOKKOS_INLINE_FUNCTION
#define KOKKOS_INLINE_FUNCTION inline
#endif

// Provide fallbacks for Kokkos math functions when not using Kokkos
#ifndef KOKKOS_CORE_HPP
namespace Kokkos {
    inline double fabs(double x) { return std::fabs(x); }
    inline double log(double x) { return std::log(x); }
    inline double pow(double x, double y) { return std::pow(x, y); }
}
#endif

namespace mpfun {
namespace detail {

// Fundamental Constants (matching MPFUN's mpfuna.f90)
constexpr int BITS_PER_WORD = 60;
constexpr int64_t RADIX = 1LL << BITS_PER_WORD;
constexpr int64_t MASK = RADIX - 1;
constexpr double DIGITS_PER_WORD = 18.061799739838871713;
constexpr double LOG_RADIX = 41.5888308335967186;
constexpr int METADATA_WORDS = 4;
constexpr int MIN_PRECISION_WORDS = 4;
constexpr int64_t MAX_EXPONENT = (1LL << 31) / BITS_PER_WORD;
constexpr int HALF_BITS = BITS_PER_WORD / 2;
constexpr int64_t HALF_MASK = (1LL << HALF_BITS) - 1;
constexpr int64_t HALF_RADIX = 1LL << HALF_BITS;
constexpr double COMPARE_FUZZ = 1.0 / (1LL << 50);
constexpr double LOG2_E = 1.4426950408889633;

// Array Index Constants
constexpr int IDX_ALLOCATED = 0;
constexpr int IDX_PRECISION = 1;
constexpr int IDX_SIGN_LENGTH = 2;
constexpr int IDX_EXPONENT = 3;
constexpr int IDX_MANTISSA_START = 4;

// Helper Functions
KOKKOS_INLINE_FUNCTION
int extract_sign(int64_t sign_length) {
    if (sign_length > 0) return 1;
    if (sign_length < 0) return -1;
    return 0;
}

KOKKOS_INLINE_FUNCTION
int extract_length(int64_t sign_length) {
    return static_cast<int>(sign_length > 0 ? sign_length : -sign_length);
}

KOKKOS_INLINE_FUNCTION
int64_t combine_sign_length(int sign, int length) {
    return static_cast<int64_t>(sign) * static_cast<int64_t>(length);
}

KOKKOS_INLINE_FUNCTION
int64_t arith_shift_right(int64_t x, int n) {
    return x >> n;
}

KOKKOS_INLINE_FUNCTION
int64_t shift_left(int64_t x, int n) {
    return x << n;
}

KOKKOS_INLINE_FUNCTION
void split_word(int64_t word, int64_t& high, int64_t& low) {
    high = arith_shift_right(word, HALF_BITS);
    low = word - shift_left(high, HALF_BITS);
}

} // namespace detail
} // namespace mpfun

#endif // MPFUN_CORE_REPRESENTATION_HPP
