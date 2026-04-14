/**
 * @file config.hpp
 * @brief Build configuration for Kokkos-MPFloat
 * 
 * This header contains compile-time configuration options for the library.
 * 
 * Copyright (c) 2024 Argonne National Laboratory
 */

#ifndef MPFUN_CONFIG_HPP
#define MPFUN_CONFIG_HPP

// Include Kokkos if available, otherwise use fallback macros and math wrappers
#if __has_include(<Kokkos_Core.hpp>)
#include <Kokkos_Core.hpp>
#define MPFUN_HAS_KOKKOS 1
#else
#define MPFUN_HAS_KOKKOS 0
// Fallback KOKKOS_INLINE_FUNCTION when Kokkos is not available
#ifndef KOKKOS_INLINE_FUNCTION
#define KOKKOS_INLINE_FUNCTION inline
#endif
#endif

// Provide Kokkos math function wrappers when not using Kokkos
// These are used by transcendental functions
#if !MPFUN_HAS_KOKKOS
#include <cmath>
namespace Kokkos {
    inline double fabs(double x) { return std::fabs(x); }
    inline double log(double x) { return std::log(x); }
    inline double pow(double x, double y) { return std::pow(x, y); }
    inline double sqrt(double x) { return std::sqrt(x); }
    inline double atan2(double y, double x) { return std::atan2(y, x); }
    inline double sin(double x) { return std::sin(x); }
    inline double cos(double x) { return std::cos(x); }
}
#endif

namespace mpfun {

/// Library version
constexpr int VERSION_MAJOR = 0;
constexpr int VERSION_MINOR = 1;
constexpr int VERSION_PATCH = 0;

/// Default precision in mantissa words (~100 digits)
constexpr int DEFAULT_PRECISION_WORDS = 6;

/// Maximum supported precision in words
constexpr int MAX_PRECISION_WORDS = 200000;

/// Whether to enable debug assertions
#ifdef NDEBUG
constexpr bool DEBUG_MODE = false;
#else
constexpr bool DEBUG_MODE = true;
#endif

/// FFT crossover point: use FFT-based multiplication above this many words
/// Default 1111 words corresponds to approximately 20,000 digits
constexpr int FFT_CROSSOVER_WORDS = 1111;

} // namespace mpfun

#endif // MPFUN_CONFIG_HPP
