# Transcendental Function Algorithm Analysis

**Author:** Agent 2 (Transcendentals Lead)  
**Date:** 2026-04-14  
**Source:** MPFUN2020-Fort v33, `mpfund.f90`

This document analyzes the algorithms used in MPFUN2020 for transcendental functions, documenting their implementation strategy, dependencies, and convergence methods.

---

## Table of Contents

1. [Overview](#overview)
2. [Constants Required](#constants-required)
3. [Function Analysis](#function-analysis)
   - [mpexp - Exponential](#mpexp---exponential)
   - [mplog - Natural Logarithm](#mplog---natural-logarithm)
   - [mpcssn - Sine and Cosine](#mpcssn---sine-and-cosine)
   - [mpang - Arctangent/Angle Computation](#mpang---arctangentangle-computation)
   - [mpcssh - Hyperbolic Sine and Cosine](#mpcssh---hyperbolic-sine-and-cosine)
4. [Supporting Functions](#supporting-functions)
5. [Dependency Graph](#dependency-graph)
6. [Precision Management](#precision-management)
7. [Port Considerations](#port-considerations)

---

## Overview

MPFUN2020's transcendental functions employ several key strategies:

1. **Argument reduction**: Reduce arguments to small ranges where Taylor series converge quickly
2. **Dynamic precision**: Adjust working precision during Newton iterations
3. **Series acceleration**: Divide reduced arguments by powers of 2, then use double-angle formulas
4. **Newton-Raphson for inverse functions**: Solve f(x) = a using Newton's method when Taylor series converge slowly

All functions depend on precomputed constants (π, log(2), Euler's gamma) stored in module arrays.

---

## Constants Required

| Constant | Module Array | Usage |
|----------|--------------|-------|
| π | `mppicon` | Trigonometric reduction, angle computation |
| log(2) | `mplog2con` | Exponential argument reduction |
| Euler's γ | `mpegammacon` | Special functions (optional for transcendentals) |

### Constant Computation

- **π (`mppiq`)**: Salamin-Brent algorithm using AGM iterations
  ```
  A_0 = 1, B_0 = 1/√2, D_0 = √2 - 1/2
  A_k = (A_{k-1} + B_{k-1})/2
  B_k = √(A_{k-1} * B_{k-1})
  D_k = D_{k-1} - 2^k * (A_k - B_k)²
  π = (A_k + B_k)² / D_k
  ```
  - Quadratic convergence (doubles digits per iteration)
  - Core ops: add, mul, sqrt, div

- **log(2) (`mplog2q`)**: AGM-based algorithm
  ```
  log(2) = π / (2 * AGM(1, 4/2^(n/2)))
  ```
  - Where n = bits of precision
  - Requires π precomputed first
  - Core ops: AGM (add, mul, sqrt), div

---

## Function Analysis

### mpexp - Exponential

**Location:** `mpfund.f90:mpexp` (lines ~715-820)

#### Algorithm

1. **Overflow/Underflow check**: Validate input magnitude
2. **Argument reduction**: 
   ```
   t = a - log(2) * nint(a / log(2))
   nz = nint(a / log(2))
   ```
   This reduces t to [-log(2)/2, log(2)/2]

3. **Further reduction**: Divide by 2^NQ where NQ = (mpnw * mpnbt)^0.4
   ```
   s = t / 2^NQ
   ```

4. **Taylor series** for exp(s):
   ```
   exp(s) = 1 + s + s²/2! + s³/3! + s⁴/4! + ...
   ```
   
5. **Recovery**: Square the result NQ times
   ```
   exp(t) = exp(s)^(2^NQ)
   ```

6. **Final adjustment**: Multiply by 2^NZ
   ```
   exp(a) = exp(t) * 2^NZ
   ```

#### Core Operations Required

| Operation | Function | Notes |
|-----------|----------|-------|
| Division | `mpdiv` | a / log(2) |
| Nearest integer | `mpnint` | Round to nearest |
| Multiplication | `mpmul`, `mpmuld` | Taylor terms, squaring |
| Addition | `mpadd` | Series accumulation |
| Scalar division | `mpdivd` | t / 2^NQ, term / k |

#### Constants Needed

- `log(2)` - precomputed in `mplog2con`

#### Convergence Strategy

- Taylor series with **linear precision reduction**: As terms get smaller, working precision decreases
- Condition: `s2(2) == 0 .or. s2(3) < s0(3) - mpnw1` (term is zero or negligible)
- Iteration limit: 1,000,000 (safety bound)

#### Notes for Port

- NQ scaling factor ~30-50 for typical precision levels (1000-10000 digits)
- The squaring loop is the most expensive part for high precision
- Consider fused multiply-add opportunities

---

### mplog - Natural Logarithm

**Location:** `mpfund.f90:mplog` (lines ~825-940)

#### Algorithm

Unlike exp, **log does NOT use Taylor series** directly (too slow convergence). Instead:

1. **Special case**: If input ≈ 1, use Taylor series for log(1+x)
   ```
   log(1+x) = x - x²/2 + x³/3 - x⁴/4 + ...
   ```
   Used when |a - 1| is very small (exponent < -2)

2. **General case**: Newton-Raphson iteration to solve exp(b) = a
   ```
   x_{k+1} = x_k + [a - exp(x_k)] / exp(x_k)
   ```
   
3. **Initial approximation**: Double-precision log
   ```
   t1 = log(mantissa) + exponent * log(2)
   ```

4. **Dynamic precision**: Start at mpnw1=4 words, double each iteration
   - Number of iterations: MQ = ceil(log₂(mpnw))

#### Core Operations Required

| Operation | Function | Notes |
|-----------|----------|-------|
| Exponential | `mpexp` | Core of Newton iteration |
| Subtraction | `mpsub` | a - exp(x_k) |
| Division | `mpdiv` | Newton update |
| Addition | `mpadd` | x + correction |
| Multiplication | `mpmul` | Taylor series terms |

#### Constants Needed

- `log(2)` - for initial approximation
- None directly, but `mpexp` needs `log(2)`

#### Convergence Strategy

- **Quadratic convergence** for Newton-Raphson (digits double each iteration)
- NIT = 3 extra iterations after reaching full precision for stability
- Special Taylor series path for arguments very close to 1

#### Critical Insight

Log is expensive because it calls exp internally. For the port:
- Optimize exp first
- Consider caching intermediate exp values if computing many logs

---

### mpcssn - Sine and Cosine

**Location:** `mpfund.f90:mpcssnr` (lines ~430-580)

This computes **both sin and cos together** (more efficient than separate).

#### Algorithm

1. **Zero check**: If input is 0, return cos=1, sin=0

2. **Verify constants**: Check π is precomputed

3. **Argument reduction to [-π, π]**:
   ```
   t = a - 2π * nint(a / 2π)
   ```

4. **Further reduction**: Divide by 2^NQ where NQ = sqrt(0.5 * mpnw * mpnbt)
   ```
   s = t / 2^NQ
   ```
   Skip if argument already very small (exponent < -1)

5. **Taylor series for sin(s)**:
   ```
   sin(s) = s - s³/3! + s⁵/5! - s⁷/7! + ...
   ```

6. **Recovery using double-angle formulas**:
   First iteration (we have sin, need cos):
   ```
   cos(2x) = 1 - 2*sin²(x)   [computed as 2*(1/2 - sin²(x))]
   ```
   
   Subsequent iterations:
   ```
   cos(2x) = 2*cos²(x) - 1   [computed as 2*(cos²(x) - 1/2)]
   ```

7. **Compute sin from cos**:
   ```
   sin(t) = ±√(1 - cos²(t))
   ```
   Sign determined by sign of original reduced sin

#### Core Operations Required

| Operation | Function | Notes |
|-----------|----------|-------|
| Division | `mpdiv`, `mpdivd` | Argument reduction, Taylor terms |
| Multiplication | `mpmul`, `mpmuld` | Taylor series, double-angle |
| Subtraction | `mpsub` | Taylor series, 1 - cos² |
| Addition | `mpadd` | Series accumulation |
| Square root | `mpsqrt` | Final sin computation |
| Nearest integer | `mpnint` | Argument reduction |

#### Constants Needed

- `π` (precomputed in `mppicon`)
- Literal constants: 1, 1/2 (constructed inline)

#### Convergence Strategy

- Taylor series with linear precision reduction
- NQ chosen to balance: more divisions vs more double-angle iterations
- Typical NQ ≈ 50-100 for high precision

#### Design Decision: Joint Computation

MPFUN computes sin and cos together because:
1. They share the reduced argument computation
2. sin is derived from cos via sqrt (cheaper than separate Taylor)
3. Many applications need both (e.g., complex exp)

**Recommendation:** Keep this design in the port as `sincos()`.

---

### mpang - Arctangent/Angle Computation

**Location:** `mpfund.f90:mpang` (lines ~140-280)

This computes the angle (atan2 equivalent) for a point (x, y) in the plane.

#### Algorithm

1. **Special cases**:
   - Both x, y = 0: Error
   - x = 0: Return ±π/2
   - y = 0: Return 0 or π

2. **Normalize**: Compute x' = x/r, y' = y/r where r = √(x² + y²)
   This ensures x'² + y'² = 1

3. **Initial approximation**: Double-precision atan2(y', x')

4. **Select Newton iteration** based on which is smaller:
   - If |x| ≤ |y|: Solve cos(a) = x using
     ```
     z_{k+1} = z_k - [x - cos(z_k)] / sin(z_k)
     ```
   - If |y| < |x|: Solve sin(a) = y using
     ```
     z_{k+1} = z_k + [y - sin(z_k)] / cos(z_k)
     ```
   
   Choosing the larger denominator improves numerical stability.

5. **Dynamic precision**: Start small, double each iteration

#### Core Operations Required

| Operation | Function | Notes |
|-----------|----------|-------|
| Sin/Cos | `mpcssnr` | Newton iteration |
| Division | `mpdiv` | Normalization, Newton update |
| Multiplication | `mpmul` | x², y², normalization |
| Addition/Subtraction | `mpadd`, `mpsub` | x² + y², corrections |
| Square root | `mpsqrt` | Normalization |
| Scalar multiply | `mpmuld` | π/2, -π/2 |

#### Constants Needed

- `π` (for special cases and quadrant handling)

#### Convergence Strategy

- Newton-Raphson with quadratic convergence
- MQ = ceil(log₂(mpnw)) iterations
- NIT = 3 extra iterations for robustness

#### Implementation Note

This is essentially `atan2(y, x)` with proper quadrant handling built in. The port should expose both:
- `atan(x)` → call `mpang(1, x)` 
- `atan2(y, x)` → call `mpang(x, y)`

---

### mpcssh - Hyperbolic Sine and Cosine

**Location:** `mpfund.f90:mpcsshr` (lines ~350-430)

Computes both sinh and cosh together.

#### Algorithm

Two paths based on argument magnitude:

**Path 1: Small argument (exponent < -1)**

Use Taylor series for sinh (avoids precision loss from exp):
```
sinh(s) = s + s³/3! + s⁵/5! + s⁷/7! + ...
```

Then compute cosh from sinh:
```
cosh(x) = √(1 + sinh²(x))
```

**Path 2: General case**

Use exponential definition:
```
t = exp(a)
cosh(a) = (t + 1/t) / 2
sinh(a) = (t - 1/t) / 2
```

#### Core Operations Required

| Operation | Function | Notes |
|-----------|----------|-------|
| Exponential | `mpexp` | General case |
| Division | `mpdiv`, `mpdivd` | 1/t, Taylor terms |
| Multiplication | `mpmul`, `mpmuld` | Taylor series |
| Addition/Subtraction | `mpadd`, `mpsub` | Exp formula |
| Square root | `mpsqrt` | Small arg: cosh from sinh |

#### Constants Needed

- None directly (exp needs log(2))

#### Convergence Strategy

- Small argument: Taylor series with precision reduction
- General case: Single exp call, no iteration

#### Design Note

Like sin/cos, computing both together is more efficient. The port should provide:
- `sinhcosh(x, &sinh_out, &cosh_out)` - joint computation
- `sinh(x)`, `cosh(x)` - convenience wrappers

---

## Supporting Functions

### mpagmr - Arithmetic-Geometric Mean

**Location:** `mpfund.f90:mpagmr` (lines ~50-100)

Used for computing π and log(2).

```
a_0 = a, b_0 = b
a_{k+1} = (a_k + b_k) / 2
b_{k+1} = √(a_k * b_k)
```

Iterate until a_k = b_k to working precision.

**Core ops:** add, mul, sqrt, div by 2

### mppiq - Pi Computation

Salamin-Brent quadratically convergent algorithm. See [Constants Required](#constants-required).

### mplog2q - Log(2) Computation

AGM-based. Requires π first.

---

## Dependency Graph

```
                    ┌─────────────┐
                    │   mppiq     │ (π)
                    │ (AGM-based) │
                    └──────┬──────┘
                           │
              ┌────────────┴────────────┐
              ▼                         ▼
       ┌─────────────┐           ┌─────────────┐
       │  mplog2q    │           │  mpcssnr    │
       │  (log 2)    │           │ (sin/cos)   │
       └──────┬──────┘           └──────┬──────┘
              │                         │
              ▼                         ▼
       ┌─────────────┐           ┌─────────────┐
       │   mpexp     │◄──────────│   mpang     │
       │    (exp)    │           │  (atan2)    │
       └──────┬──────┘           └─────────────┘
              │                         ▲
              ▼                         │
       ┌─────────────┐                  │
       │   mplog     │──────────────────┘
       │    (log)    │      (mplog uses mpexp)
       └──────┬──────┘
              │
              ▼
       ┌─────────────┐
       │   mpcssh    │
       │ (sinh/cosh) │
       └──────┬──────┘
              │
              ▼
       ┌─────────────┐
       │  mppower    │
       │   (a^b)     │
       └─────────────┘

Core arithmetic dependencies (all transcendentals):
┌─────────────────────────────────────────────┐
│  mpadd, mpsub, mpmul, mpdiv, mpsqrt        │
│  mpmuld, mpdivd (scalar operations)         │
│  mpnint, mpcpr, mpeq (utilities)            │
└─────────────────────────────────────────────┘
```

---

## Precision Management

### Dynamic Precision Pattern

Most Newton-Raphson iterations use this pattern:

```fortran
mpnw1 = 4                    ! Start with minimal precision
do k = 1, mq
  if (k > 1) mpnw1 = min(2 * mpnw1 - 2, mpnw) + 1   ! Double precision
  
  ! ... computation at precision mpnw1 ...
  
  if (k == mq - nit .and. iq == 0) then   ! Extra iterations at full precision
    iq = 1
    goto repeat_iteration
  endif
enddo
```

### Linear Precision Reduction for Series

```fortran
mpnw2 = mpnw1
do j = 1, itrmx
  ! ... compute next term at precision mpnw2 ...
  
  ! Check convergence
  if (term_is_negligible) exit
  
  ! Reduce precision for next term
  mpnw2 = min(max(mpnw1 + int(term_exp - sum_exp) + 1, 4), mpnw1)
enddo
```

This optimization reduces computation as terms become negligible.

---

## Port Considerations

### Blocking on Agent 1

The following core operations must be available before implementing transcendentals:
- `mpadd`, `mpsub` - Addition/subtraction
- `mpmul`, `mpmuld` - Multiplication (full and scalar)
- `mpdiv`, `mpdivd` - Division (full and scalar)
- `mpsqrt` - Square root
- `mpcpr` - Comparison
- `mpnint` - Nearest integer
- `mpeq` - Copy/assignment

### Recommended Implementation Order

1. **Constants first**: π, log(2) (required by most functions)
2. **exp** - Foundation, doesn't depend on other transcendentals
3. **log** - Depends on exp
4. **sincos** - Independent path from exp/log
5. **atan2/atan** - Depends on sincos
6. **sinhcosh** - Depends on exp
7. **Derived functions**: tan, tanh, asin, acos, asinh, acosh, atanh

### KOKKOS_INLINE_FUNCTION Considerations

1. **No dynamic allocation**: All temporaries use fixed-size arrays
2. **No exceptions**: Use NaN or error codes for invalid input
3. **Recursion**: Avoid deep recursion (Newton iterations are loops, not recursive)
4. **Loop bounds**: All iteration limits are explicit constants

### Potential Optimizations

1. **Fused operations**: Some algorithms could benefit from FMA
2. **Precomputation**: π, log(2) can be precomputed at compile time for common precisions
3. **Parallelism**: Individual function calls are serial, but many independent calls can be batched
4. **Precision hints**: If user knows result precision needed, could skip extra iterations

---

## References

1. Bailey, D.H. "MPFUN2020: A new thread-safe arbitrary precision package"
2. Brent, R.P. "Fast Multiple-Precision Evaluation of Elementary Functions"
3. Salamin, E. "Computation of π Using Arithmetic-Geometric Mean"
