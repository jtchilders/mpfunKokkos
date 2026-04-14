#ifndef MPFUN_FWD_HPP
#define MPFUN_FWD_HPP

namespace mpfun {

// Forward declaration of the main type
template <int WORDS>
class MPFloat;

// Common precision level aliases (forward declared)
// Actual definitions in mp_float.hpp
using MPFloat50   = MPFloat<3>;    // ~50 digits
using MPFloat100  = MPFloat<6>;    // ~100 digits
using MPFloat200  = MPFloat<12>;   // ~200 digits
using MPFloat500  = MPFloat<28>;   // ~500 digits
using MPFloat1000 = MPFloat<56>;   // ~1000 digits

// Short aliases
using mp50   = MPFloat50;
using mp100  = MPFloat100;
using mp200  = MPFloat200;
using mp500  = MPFloat500;
using mp1000 = MPFloat1000;

} // namespace mpfun

#endif // MPFUN_FWD_HPP
