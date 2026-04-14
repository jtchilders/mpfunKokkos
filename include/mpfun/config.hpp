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

// Include Kokkos if available, otherwise use fallback macros
#if __has_include(<Kokkos_Core.hpp>)
#include <Kokkos_Core.hpp>
#else
// Fallback when Kokkos is not available
#ifndef KOKKOS_INLINE_FUNCTION
#define KOKKOS_INLINE_FUNCTION inline
#endif
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
