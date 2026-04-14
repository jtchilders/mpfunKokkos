# Interface Contracts

This document defines the agreed-upon interfaces between agents. **Do not deviate without coordinator approval.**

---

## Memory Layout

Following MPFUN20-Fort conventions:

```
data_[0] = allocated_words     (capacity of mantissa, not counting metadata)
data_[1] = working_precision   (current precision in words)
data_[2] = sign * length       (sign: ±1, length: number of mantissa words used)
data_[3] = exponent            (power of RADIX = 2^60)
data_[4] = mantissa[0]         (most significant word)
data_[5] = mantissa[1]
...
data_[3+length] = mantissa[length-1]  (least significant word)
```

### Constants

```cpp
namespace mpfun {
namespace detail {

constexpr int BITS_PER_WORD = 60;
constexpr int64_t RADIX = 1LL << 60;           // 2^60
constexpr int64_t MASK = RADIX - 1;            // 0x0FFFFFFFFFFFFFFF
constexpr double DIGITS_PER_WORD = 18.0617997398;  // log10(2^60)
constexpr int METADATA_WORDS = 4;

} // namespace detail
} // namespace mpfun
```

### Zero Representation

A value is zero when `data_[2] == 0`. The exponent and mantissa are undefined when zero.

### Sign Convention

- Positive: `data_[2] > 0`
- Negative: `data_[2] < 0`
- Zero: `data_[2] == 0`

The magnitude (number of mantissa words) is `abs(data_[2])`.

---

## MPFloat Class Template

```cpp
namespace mpfun {

template <int WORDS>
class MPFloat {
public:
    // Type aliases
    using word_type = int64_t;
    
    // Constants
    static constexpr int capacity() { return WORDS; }
    static constexpr int approx_digits() { 
        return static_cast<int>(WORDS * detail::DIGITS_PER_WORD); 
    }
    
    // === Constructors ===
    
    KOKKOS_INLINE_FUNCTION 
    MPFloat();  // Zero-initialized
    
    KOKKOS_INLINE_FUNCTION 
    explicit MPFloat(double d);
    
    KOKKOS_INLINE_FUNCTION 
    explicit MPFloat(int64_t i);
    
    KOKKOS_INLINE_FUNCTION 
    explicit MPFloat(int i);
    
    // String constructor (host only)
    explicit MPFloat(const char* str);
    
    // === Assignment ===
    
    KOKKOS_INLINE_FUNCTION 
    MPFloat& operator=(const MPFloat& rhs) = default;
    
    KOKKOS_INLINE_FUNCTION 
    MPFloat& operator=(double d);
    
    // === Arithmetic Operators ===
    
    KOKKOS_INLINE_FUNCTION 
    MPFloat operator+(const MPFloat& rhs) const;
    
    KOKKOS_INLINE_FUNCTION 
    MPFloat operator-(const MPFloat& rhs) const;
    
    KOKKOS_INLINE_FUNCTION 
    MPFloat operator*(const MPFloat& rhs) const;
    
    KOKKOS_INLINE_FUNCTION 
    MPFloat operator/(const MPFloat& rhs) const;
    
    KOKKOS_INLINE_FUNCTION 
    MPFloat operator-() const;  // Unary negation
    
    // === Compound Assignment ===
    
    KOKKOS_INLINE_FUNCTION 
    MPFloat& operator+=(const MPFloat& rhs);
    
    KOKKOS_INLINE_FUNCTION 
    MPFloat& operator-=(const MPFloat& rhs);
    
    KOKKOS_INLINE_FUNCTION 
    MPFloat& operator*=(const MPFloat& rhs);
    
    KOKKOS_INLINE_FUNCTION 
    MPFloat& operator/=(const MPFloat& rhs);
    
    // === Comparison Operators ===
    
    KOKKOS_INLINE_FUNCTION 
    bool operator==(const MPFloat& rhs) const;
    
    KOKKOS_INLINE_FUNCTION 
    bool operator!=(const MPFloat& rhs) const;
    
    KOKKOS_INLINE_FUNCTION 
    bool operator<(const MPFloat& rhs) const;
    
    KOKKOS_INLINE_FUNCTION 
    bool operator<=(const MPFloat& rhs) const;
    
    KOKKOS_INLINE_FUNCTION 
    bool operator>(const MPFloat& rhs) const;
    
    KOKKOS_INLINE_FUNCTION 
    bool operator>=(const MPFloat& rhs) const;
    
    // === Accessors ===
    
    KOKKOS_INLINE_FUNCTION 
    int sign() const;  // Returns -1, 0, or +1
    
    KOKKOS_INLINE_FUNCTION 
    int64_t exponent() const;
    
    KOKKOS_INLINE_FUNCTION 
    int mantissa_length() const;  // Number of words used
    
    KOKKOS_INLINE_FUNCTION 
    int precision_words() const;  // Working precision
    
    KOKKOS_INLINE_FUNCTION 
    bool is_zero() const;
    
    KOKKOS_INLINE_FUNCTION 
    bool is_negative() const;
    
    KOKKOS_INLINE_FUNCTION 
    bool is_positive() const;
    
    // === Conversion ===
    
    KOKKOS_INLINE_FUNCTION 
    double to_double() const;
    
    // String conversion (host only)
    std::string to_string(int digits = -1) const;
    
    // === Raw Access (use with caution) ===
    
    KOKKOS_INLINE_FUNCTION 
    word_type* data() { return data_; }
    
    KOKKOS_INLINE_FUNCTION 
    const word_type* data() const { return data_; }
    
private:
    word_type data_[WORDS + detail::METADATA_WORDS];
    
    // Internal helpers
    KOKKOS_INLINE_FUNCTION void normalize();
    KOKKOS_INLINE_FUNCTION void set_zero();
};

} // namespace mpfun
```

---

## Free Functions

### Core Math (Agent 1)

```cpp
namespace mpfun {

template <int W>
KOKKOS_INLINE_FUNCTION 
MPFloat<W> abs(const MPFloat<W>& x);

template <int W>
KOKKOS_INLINE_FUNCTION 
MPFloat<W> sqrt(const MPFloat<W>& x);

template <int W>
KOKKOS_INLINE_FUNCTION 
MPFloat<W> nroot(const MPFloat<W>& x, int n);  // n-th root

template <int W>
KOKKOS_INLINE_FUNCTION 
int compare(const MPFloat<W>& a, const MPFloat<W>& b);  // Returns -1, 0, +1

template <int W>
KOKKOS_INLINE_FUNCTION 
MPFloat<W> max(const MPFloat<W>& a, const MPFloat<W>& b);

template <int W>
KOKKOS_INLINE_FUNCTION 
MPFloat<W> min(const MPFloat<W>& a, const MPFloat<W>& b);

} // namespace mpfun
```

### Transcendentals (Agent 2)

```cpp
namespace mpfun {

// Exponential / Logarithm
template <int W>
KOKKOS_INLINE_FUNCTION 
MPFloat<W> exp(const MPFloat<W>& x);

template <int W>
KOKKOS_INLINE_FUNCTION 
MPFloat<W> expm1(const MPFloat<W>& x);  // exp(x) - 1, accurate for small x

template <int W>
KOKKOS_INLINE_FUNCTION 
MPFloat<W> log(const MPFloat<W>& x);

template <int W>
KOKKOS_INLINE_FUNCTION 
MPFloat<W> log1p(const MPFloat<W>& x);  // log(1 + x), accurate for small x

template <int W>
KOKKOS_INLINE_FUNCTION 
MPFloat<W> log10(const MPFloat<W>& x);

template <int W>
KOKKOS_INLINE_FUNCTION 
MPFloat<W> log2(const MPFloat<W>& x);

template <int W>
KOKKOS_INLINE_FUNCTION 
MPFloat<W> pow(const MPFloat<W>& base, const MPFloat<W>& exp);

template <int W>
KOKKOS_INLINE_FUNCTION 
MPFloat<W> pow(const MPFloat<W>& base, int exp);  // Integer exponent (faster)

// Trigonometric
template <int W>
KOKKOS_INLINE_FUNCTION 
MPFloat<W> sin(const MPFloat<W>& x);

template <int W>
KOKKOS_INLINE_FUNCTION 
MPFloat<W> cos(const MPFloat<W>& x);

template <int W>
KOKKOS_INLINE_FUNCTION 
void sincos(const MPFloat<W>& x, MPFloat<W>& sin_out, MPFloat<W>& cos_out);

template <int W>
KOKKOS_INLINE_FUNCTION 
MPFloat<W> tan(const MPFloat<W>& x);

template <int W>
KOKKOS_INLINE_FUNCTION 
MPFloat<W> asin(const MPFloat<W>& x);

template <int W>
KOKKOS_INLINE_FUNCTION 
MPFloat<W> acos(const MPFloat<W>& x);

template <int W>
KOKKOS_INLINE_FUNCTION 
MPFloat<W> atan(const MPFloat<W>& x);

template <int W>
KOKKOS_INLINE_FUNCTION 
MPFloat<W> atan2(const MPFloat<W>& y, const MPFloat<W>& x);

// Hyperbolic
template <int W>
KOKKOS_INLINE_FUNCTION 
MPFloat<W> sinh(const MPFloat<W>& x);

template <int W>
KOKKOS_INLINE_FUNCTION 
MPFloat<W> cosh(const MPFloat<W>& x);

template <int W>
KOKKOS_INLINE_FUNCTION 
void sinhcosh(const MPFloat<W>& x, MPFloat<W>& sinh_out, MPFloat<W>& cosh_out);

template <int W>
KOKKOS_INLINE_FUNCTION 
MPFloat<W> tanh(const MPFloat<W>& x);

template <int W>
KOKKOS_INLINE_FUNCTION 
MPFloat<W> asinh(const MPFloat<W>& x);

template <int W>
KOKKOS_INLINE_FUNCTION 
MPFloat<W> acosh(const MPFloat<W>& x);

template <int W>
KOKKOS_INLINE_FUNCTION 
MPFloat<W> atanh(const MPFloat<W>& x);

} // namespace mpfun
```

### Constants (shared responsibility)

```cpp
namespace mpfun {

// Returns pi to the precision of MPFloat<W>
template <int W>
KOKKOS_INLINE_FUNCTION 
MPFloat<W> pi();

// Returns e (Euler's number) to the precision of MPFloat<W>
template <int W>
KOKKOS_INLINE_FUNCTION 
MPFloat<W> e();

// Returns log(2) to the precision of MPFloat<W>
template <int W>
KOKKOS_INLINE_FUNCTION 
MPFloat<W> ln2();

// Returns Euler-Mascheroni constant
template <int W>
KOKKOS_INLINE_FUNCTION 
MPFloat<W> euler_gamma();

} // namespace mpfun
```

---

## Type Aliases (Convenience)

```cpp
namespace mpfun {

// Common precision levels
using MPFloat50   = MPFloat<3>;    // ~50 digits
using MPFloat100  = MPFloat<6>;    // ~100 digits
using MPFloat200  = MPFloat<12>;   // ~200 digits
using MPFloat500  = MPFloat<28>;   // ~500 digits
using MPFloat1000 = MPFloat<56>;   // ~1000 digits

// Shorter aliases
using mp50   = MPFloat50;
using mp100  = MPFloat100;
using mp200  = MPFloat200;
using mp500  = MPFloat500;
using mp1000 = MPFloat1000;

} // namespace mpfun
```

---

## Error Handling

Since we're targeting GPU code (`KOKKOS_INLINE_FUNCTION`), we cannot use exceptions. 

### Strategy

1. **Precondition checking**: Validate inputs in debug builds only
2. **NaN propagation**: Invalid operations produce a "NaN" state
3. **Assertions**: Use `KOKKOS_ASSERT` for debug builds

### NaN Representation

A special invalid/NaN state is indicated by:
```cpp
data_[0] = -1;  // Indicates NaN
```

Check with:
```cpp
KOKKOS_INLINE_FUNCTION bool is_nan() const { return data_[0] < 0; }
```

---

## File Organization

```
include/mpfun/
├── mpfun.hpp              # Single-include umbrella header
├── config.hpp             # Build configuration macros
├── fwd.hpp                # Forward declarations
├── mp_float.hpp           # MPFloat class definition
├── core/
│   ├── representation.hpp # Memory layout, constants
│   ├── add.hpp            # Addition implementation
│   ├── sub.hpp            # Subtraction implementation  
│   ├── mul.hpp            # Multiplication implementation
│   ├── div.hpp            # Division implementation
│   ├── sqrt.hpp           # Square root implementation
│   └── compare.hpp        # Comparison implementation
├── transcendental/
│   ├── exp.hpp
│   ├── log.hpp
│   ├── trig.hpp           # sin, cos, tan
│   ├── inverse_trig.hpp   # asin, acos, atan
│   └── hyperbolic.hpp     # sinh, cosh, tanh, etc.
├── constants/
│   ├── pi.hpp
│   ├── e.hpp
│   └── precomputed.hpp    # Tables for high precision
└── io/                    # Host-only I/O
    ├── string.hpp
    └── stream.hpp
```

---

## Testing Requirements

Each function must have:

1. **Basic correctness test**: Known input → expected output
2. **Edge cases**: Zero, negative, very large, very small
3. **Precision validation**: Compare against Fortran MPFUN output
4. **Identity tests**: e.g., `sin²(x) + cos²(x) = 1`

### Test Naming Convention

```
test_<function>_<case>.cpp
```

Examples:
- `test_add_basic.cpp`
- `test_add_overflow.cpp`
- `test_sin_special_values.cpp`

---

## Versioning

API changes require:
1. Update to this document
2. Notification to all agents
3. Coordinator approval

Current interface version: **1.0-draft**
