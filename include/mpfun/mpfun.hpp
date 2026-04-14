/**
 * @file mpfun.hpp
 * @brief Single-include umbrella header for Kokkos-MPFloat
 * 
 * This header includes all components of the Kokkos-MPFloat library.
 * For most uses, include only this header.
 * 
 * Example usage:
 * @code
 *   #include <mpfun/mpfun.hpp>
 *   
 *   using namespace mpfun;
 *   
 *   MPFloat100 a(3.14159);
 *   MPFloat100 b(2.71828);
 *   MPFloat100 c = a * b + a / b;
 * @endcode
 * 
 * Copyright (c) 2024 Argonne National Laboratory
 * Ported from MPFUN2020-Fort by David H. Bailey
 */

#ifndef MPFUN_MPFUN_HPP
#define MPFUN_MPFUN_HPP

// Configuration
#include "config.hpp"

// Core types
#include "mp_float.hpp"

// Note: mp_float.hpp includes all core operation headers:
// - core/representation.hpp
// - core/add.hpp
// - core/mul.hpp  
// - core/div.hpp
// - core/sqrt.hpp
// - core/compare.hpp

// TODO: Add when implemented
// #include "mp_complex.hpp"
// #include "transcendental/exp.hpp"
// #include "transcendental/log.hpp"
// #include "transcendental/trig.hpp"
// #include "constants/pi.hpp"

#endif // MPFUN_MPFUN_HPP
