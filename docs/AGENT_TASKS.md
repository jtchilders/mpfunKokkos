# Agent Task Assignments

## Overview

Three parallel agents plus coordinator (main session). All work in `~/projects/bailey-port/`.

---

## Agent 1: Core Arithmetic Lead

**Label**: `mpfun-core`  
**Priority**: Critical path — others depend on this  
**Estimated Duration**: 1-2 weeks

### Primary Responsibilities

1. **Define the MPFloat type** (`include/mpfun/mp_float.hpp`)
   - Template class `MPFloat<WORDS>`
   - Memory layout matching MPFUN conventions
   - All methods marked `KOKKOS_INLINE_FUNCTION`

2. **Implement core operations** (in `include/mpfun/core/`)
   - `representation.hpp` — data layout, accessors, metadata handling
   - `add_sub.hpp` — `mpadd`, `mpsub` (with normalization)
   - `mul.hpp` — `mpmul` (schoolbook first, FFT later)
   - `div.hpp` — `mpdiv`, `mpdivd`
   - `sqrt.hpp` — `mpsqrt`
   - `compare.hpp` — `mpcpr`, comparison operators

3. **Port utility functions**
   - `mpnorm` — normalize after operations
   - `mpeq` — copy/assignment
   - `mpabs` — absolute value
   - `mpneg` — negation

### Key Files to Analyze

```
mpfun20-fort-v33/fortran-var1/
├── mpfuna.f90   # Parameters, constants — understand data layout
├── mpfunb.f90   # Core arithmetic — PRIMARY PORTING TARGET
```

### Deliverables

1. Working `MPFloat<N>` with +, -, *, /, sqrt, comparisons
2. Passes basic arithmetic tests
3. Documentation of memory layout decisions

### Interface Contract (to be agreed upon)

```cpp
template <int WORDS>
class MPFloat {
public:
    static constexpr int BITS_PER_WORD = 60;
    static constexpr int64_t RADIX = 1LL << 60;
    
    KOKKOS_INLINE_FUNCTION MPFloat();
    KOKKOS_INLINE_FUNCTION explicit MPFloat(double d);
    KOKKOS_INLINE_FUNCTION explicit MPFloat(int64_t i);
    
    // Arithmetic
    KOKKOS_INLINE_FUNCTION MPFloat operator+(const MPFloat& rhs) const;
    KOKKOS_INLINE_FUNCTION MPFloat operator-(const MPFloat& rhs) const;
    KOKKOS_INLINE_FUNCTION MPFloat operator*(const MPFloat& rhs) const;
    KOKKOS_INLINE_FUNCTION MPFloat operator/(const MPFloat& rhs) const;
    KOKKOS_INLINE_FUNCTION MPFloat operator-() const;
    
    // Compound assignment
    KOKKOS_INLINE_FUNCTION MPFloat& operator+=(const MPFloat& rhs);
    KOKKOS_INLINE_FUNCTION MPFloat& operator-=(const MPFloat& rhs);
    KOKKOS_INLINE_FUNCTION MPFloat& operator*=(const MPFloat& rhs);
    KOKKOS_INLINE_FUNCTION MPFloat& operator/=(const MPFloat& rhs);
    
    // Comparison
    KOKKOS_INLINE_FUNCTION bool operator==(const MPFloat& rhs) const;
    KOKKOS_INLINE_FUNCTION bool operator!=(const MPFloat& rhs) const;
    KOKKOS_INLINE_FUNCTION bool operator<(const MPFloat& rhs) const;
    KOKKOS_INLINE_FUNCTION bool operator<=(const MPFloat& rhs) const;
    KOKKOS_INLINE_FUNCTION bool operator>(const MPFloat& rhs) const;
    KOKKOS_INLINE_FUNCTION bool operator>=(const MPFloat& rhs) const;
    
    // Accessors
    KOKKOS_INLINE_FUNCTION int sign() const;
    KOKKOS_INLINE_FUNCTION int exponent() const;
    KOKKOS_INLINE_FUNCTION int precision_words() const;
    KOKKOS_INLINE_FUNCTION bool is_zero() const;
    
    // Conversion
    KOKKOS_INLINE_FUNCTION double to_double() const;
    
private:
    int64_t data_[WORDS + 4];
    
    KOKKOS_INLINE_FUNCTION void normalize();
};

// Free functions
template <int W>
KOKKOS_INLINE_FUNCTION MPFloat<W> abs(const MPFloat<W>& x);

template <int W>
KOKKOS_INLINE_FUNCTION MPFloat<W> sqrt(const MPFloat<W>& x);
```

---

## Agent 2: Transcendentals

**Label**: `mpfun-transcendental`  
**Priority**: Blocked on Agent 1 core types  
**Estimated Duration**: 1-1.5 weeks (after Agent 1 delivers)

### Primary Responsibilities

1. **Exponential and Logarithm** (`include/mpfun/transcendental/exp_log.hpp`)
   - `exp(MPFloat)` — port `mpexp`
   - `log(MPFloat)` — port `mplog`
   - `log10(MPFloat)` — derived from log

2. **Trigonometric** (`include/mpfun/transcendental/trig.hpp`)
   - `sin(MPFloat)`, `cos(MPFloat)` — port `mpcssn` (computes both)
   - `tan(MPFloat)` — port `mptan`
   - `asin`, `acos`, `atan`, `atan2` — port `mpang` etc.

3. **Hyperbolic** (`include/mpfun/transcendental/hyperbolic.hpp`)
   - `sinh(MPFloat)`, `cosh(MPFloat)` — port `mpcssh`
   - `tanh(MPFloat)`
   - `asinh`, `acosh`, `atanh`

### Key Files to Analyze

```
mpfun20-fort-v33/fortran-var1/
├── mpfund.f90   # Transcendentals — PRIMARY PORTING TARGET
├── mpfuna.f90   # Precomputed constants (pi, log2)
```

### Initial Work (Before Agent 1 Delivers)

1. **Study the algorithms** in `mpfund.f90`
2. **Document dependencies**: What core operations does each transcendental need?
3. **Design test cases**: Known values (sin(π/6) = 0.5, exp(0) = 1, etc.)
4. **Stub out headers** with function signatures

### Deliverables

1. Working exp, log, sin, cos, tan, atan, sinh, cosh, tanh
2. All functions `KOKKOS_INLINE_FUNCTION` compatible
3. Validation tests against reference values

### Interface Contract

```cpp
namespace mpfun {

// Exponential/Logarithm
template <int W>
KOKKOS_INLINE_FUNCTION MPFloat<W> exp(const MPFloat<W>& x);

template <int W>
KOKKOS_INLINE_FUNCTION MPFloat<W> log(const MPFloat<W>& x);

template <int W>
KOKKOS_INLINE_FUNCTION MPFloat<W> log10(const MPFloat<W>& x);

template <int W>
KOKKOS_INLINE_FUNCTION MPFloat<W> pow(const MPFloat<W>& base, const MPFloat<W>& exp);

// Trigonometric
template <int W>
KOKKOS_INLINE_FUNCTION MPFloat<W> sin(const MPFloat<W>& x);

template <int W>
KOKKOS_INLINE_FUNCTION MPFloat<W> cos(const MPFloat<W>& x);

template <int W>
KOKKOS_INLINE_FUNCTION void sincos(const MPFloat<W>& x, MPFloat<W>& s, MPFloat<W>& c);

template <int W>
KOKKOS_INLINE_FUNCTION MPFloat<W> tan(const MPFloat<W>& x);

template <int W>
KOKKOS_INLINE_FUNCTION MPFloat<W> atan(const MPFloat<W>& x);

template <int W>
KOKKOS_INLINE_FUNCTION MPFloat<W> atan2(const MPFloat<W>& y, const MPFloat<W>& x);

// Hyperbolic
template <int W>
KOKKOS_INLINE_FUNCTION MPFloat<W> sinh(const MPFloat<W>& x);

template <int W>
KOKKOS_INLINE_FUNCTION MPFloat<W> cosh(const MPFloat<W>& x);

template <int W>
KOKKOS_INLINE_FUNCTION MPFloat<W> tanh(const MPFloat<W>& x);

} // namespace mpfun
```

---

## Agent 3: Test Infrastructure & Build System

**Label**: `mpfun-infra`  
**Priority**: Can start immediately, parallel to Agent 1  
**Estimated Duration**: Ongoing through project

### Primary Responsibilities

1. **Build System** (CMake)
   - CMakeLists.txt for library
   - Find/configure Kokkos
   - Optional CUDA/HIP/SYCL backends
   - Install targets

2. **Test Framework**
   - Set up GoogleTest or Catch2
   - Test harness that can run on host and device
   - Validation framework comparing against Fortran reference

3. **Reference Comparison**
   - Compile Fortran MPFUN with gfortran
   - Create test driver that outputs known values
   - Script to compare C++ output against Fortran output

4. **Benchmarking**
   - Kokkos-based benchmark harness
   - Measure: operations/sec, bandwidth, scaling

5. **CI/CD**
   - GitHub Actions workflow
   - Test matrix: CPU, CUDA (if available), multiple compilers

### Deliverables

1. Working CMake build system
2. Test harness with validation against Fortran reference
3. Benchmark framework
4. CI configuration

### Directory Structure to Create

```
kokkos-mpfun/
├── CMakeLists.txt
├── cmake/
│   ├── FindKokkos.cmake
│   └── mpfun-config.cmake.in
├── include/
│   └── mpfun/
│       └── (Agent 1, 2 headers)
├── tests/
│   ├── CMakeLists.txt
│   ├── unit/
│   │   ├── test_add.cpp
│   │   ├── test_mul.cpp
│   │   └── ...
│   ├── validation/
│   │   ├── reference_values.hpp   # From Fortran output
│   │   └── test_against_reference.cpp
│   └── fortran_reference/
│       ├── generate_reference.f90
│       └── Makefile
├── benchmarks/
│   ├── CMakeLists.txt
│   └── bench_arithmetic.cpp
├── examples/
│   ├── CMakeLists.txt
│   └── basic_usage.cpp
└── .github/
    └── workflows/
        └── ci.yml
```

### Reference Value Generation

Create Fortran program to output test values:

```fortran
! generate_reference.f90
program generate_reference
    use mpmodule
    implicit none
    type(mp_real) :: a, b, c
    
    ! Test addition
    a = '3.14159265358979323846'
    b = '2.71828182845904523536'
    c = a + b
    call mpwrite(6, 100, 80, c)
    
    ! ... more test cases
end program
```

---

## Coordinator Role (Main Session)

**You (Taylor) + Wesley in main chat**

### Responsibilities

1. **Spawn and monitor agents**
2. **Interface arbitration**: If agents disagree on API, decide
3. **Integration**: Merge work from agents, resolve conflicts
4. **Architecture decisions**: Any major changes to plan
5. **Final review**: Code review merged work

### Checkpoints

| Milestone | Expected | Verification |
|-----------|----------|--------------|
| Build system working | Day 2 | `cmake .. && make` succeeds |
| MPFloat type compiles | Day 3 | Minimal test passes |
| Add/Sub working | Day 5 | Unit tests pass |
| Mul/Div working | Day 8 | Unit tests pass |
| Phase 1 complete | Day 10 | Validation against Fortran |
| Transcendentals started | Day 8 | Agent 2 unblocked |
| Phase 2 complete | Day 15 | Full validation |

---

## Communication Guidelines

### Between Agents

- **Shared contracts**: All agents reference `docs/INTERFACES.md` (to be created)
- **File ownership**: Each agent owns specific directories, no overlap
- **Blocking issues**: Report immediately to coordinator

### Deliverable Format

Each agent delivers:
1. Working code in their assigned directories
2. Unit tests for their code
3. Documentation (comments + any necessary README)

### Merge Protocol

1. Agent commits to feature branch
2. Agent notifies coordinator: "Ready for review: [feature]"
3. Coordinator reviews, requests changes if needed
4. Coordinator merges to main

---

## Quick Reference: Fortran → C++ Translation

| Fortran | C++ Equivalent |
|---------|----------------|
| `integer(mpiknd)` | `int64_t` |
| `a(0:)` array | `data_[N]` fixed array |
| `intent(in)` | `const&` |
| `intent(out)` | return value or `&` param |
| `call mpadd(a,b,c,mpnw)` | `c = a + b` (operator overload) |
| `if (a(2) == 0)` (check zero) | `if (a.is_zero())` |
| `mpbdx = 2^60` | `static constexpr RADIX = 1LL << 60` |
