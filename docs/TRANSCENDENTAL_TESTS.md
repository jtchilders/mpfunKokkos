# Transcendental Function Test Cases

**Author:** Agent 2 (Transcendentals Lead)  
**Date:** 2026-04-14  
**Purpose:** Validation test cases for transcendental function implementations

This document specifies test cases for validating the correctness of transcendental functions. Tests are organized by function and category.

---

## Table of Contents

1. [Test Infrastructure Requirements](#test-infrastructure-requirements)
2. [Exponential and Logarithm Tests](#exponential-and-logarithm-tests)
3. [Trigonometric Tests](#trigonometric-tests)
4. [Inverse Trigonometric Tests](#inverse-trigonometric-tests)
5. [Hyperbolic Tests](#hyperbolic-tests)
6. [Inverse Hyperbolic Tests](#inverse-hyperbolic-tests)
7. [Identity Tests](#identity-tests)
8. [Edge Cases and Error Handling](#edge-cases-and-error-handling)
9. [Precision Validation](#precision-validation)
10. [Cross-Validation Against Fortran MPFUN](#cross-validation-against-fortran-mpfun)

---

## Test Infrastructure Requirements

### Tolerance Specification

For precision level W words (~18*W digits):
- **Relative tolerance**: `10^(-(18*W - 5))` — allow 5 digits of error margin
- **Absolute tolerance for near-zero**: `10^(-18*W)` — for results that should be exactly zero

### Test Macros (to be implemented by Agent 3)

```cpp
// Check that |computed - expected| / |expected| < tol
MPFUN_ASSERT_NEAR_REL(computed, expected, rel_tol);

// Check that |computed - expected| < tol  
MPFUN_ASSERT_NEAR_ABS(computed, expected, abs_tol);

// Check that computed represents NaN/invalid
MPFUN_ASSERT_NAN(computed);

// Check exact equality (for special values like 0, 1)
MPFUN_ASSERT_EQ(computed, expected);
```

---

## Exponential and Logarithm Tests

### exp(x) - Basic Values

| Input | Expected Output | Notes |
|-------|-----------------|-------|
| `0` | `1` | exp(0) = 1 exactly |
| `1` | `e = 2.71828182845904523536...` | Euler's number |
| `-1` | `1/e = 0.36787944117144232159...` | |
| `0.5` | `sqrt(e) = 1.64872127070012814684...` | |
| `log(2)` | `2` | exp(log(2)) = 2 exactly |
| `log(10)` | `10` | |
| `10` | `22026.46579480671651695...` | Large value |
| `-10` | `0.00004539992976248485154...` | Small value |
| `100` | `2.688117141816135...e43` | Very large |
| `-100` | `3.720075976020836...e-44` | Very small |

### exp(x) - Edge Cases

| Input | Expected | Notes |
|-------|----------|-------|
| `+∞` (overflow input) | Error/NaN | Input too large |
| `-∞` (underflow) | `0` | Graceful underflow to zero |
| Very small positive | `1 + x + O(x²)` | Tests expm1 accuracy |

### log(x) - Basic Values

| Input | Expected Output | Notes |
|-------|-----------------|-------|
| `1` | `0` | log(1) = 0 exactly |
| `e` | `1` | By definition |
| `e²` | `2` | |
| `2` | `0.69314718055994530941...` | log(2) |
| `10` | `2.30258509299404568401...` | log(10) |
| `0.5` | `-0.69314718055994530941...` | -log(2) |
| `1/e` | `-1` | |
| `exp(1)` | `1` | log(exp(1)) = 1 |

### log(x) - Edge Cases

| Input | Expected | Notes |
|-------|----------|-------|
| `0` | Error/NaN | log(0) undefined |
| `-1` | Error/NaN | log of negative undefined |
| `1 + 1e-50` | `~1e-50` | Tests log1p accuracy |
| Very large | `log(x)` | Tests argument reduction |

### expm1(x) and log1p(x) - Accuracy Tests

| Function | Input | Expected | Notes |
|----------|-------|----------|-------|
| `expm1` | `1e-20` | `1e-20 + 5e-41 + ...` | Near machine epsilon |
| `expm1` | `1e-50` | `1e-50 + ...` | Very small |
| `log1p` | `1e-20` | `1e-20 - 5e-41 + ...` | |
| `log1p` | `1e-50` | `1e-50 - ...` | |

---

## Trigonometric Tests

### sin(x) and cos(x) - Special Values

| x | sin(x) | cos(x) | Notes |
|---|--------|--------|-------|
| `0` | `0` | `1` | Exact |
| `π/6` | `0.5` | `sqrt(3)/2 = 0.866025...` | 30° |
| `π/4` | `sqrt(2)/2 = 0.707106...` | `sqrt(2)/2` | 45° |
| `π/3` | `sqrt(3)/2` | `0.5` | 60° |
| `π/2` | `1` | `0` | 90°, exact |
| `π` | `0` | `-1` | 180°, exact |
| `3π/2` | `-1` | `0` | 270° |
| `2π` | `0` | `1` | 360° |
| `-π/2` | `-1` | `0` | |

### sin(x) and cos(x) - Large Arguments

| x | Expected | Notes |
|---|----------|-------|
| `1e6` | Compute via reduction | Tests argument reduction |
| `1e12` | Compute via reduction | Stress test |
| `1e50` | May lose precision | Document precision limits |

### tan(x) - Values

| x | tan(x) | Notes |
|---|--------|-------|
| `0` | `0` | |
| `π/4` | `1` | Exact |
| `π/3` | `sqrt(3) = 1.732050...` | |
| `π/6` | `1/sqrt(3) = 0.577350...` | |
| `π/2` | Error/NaN or ±∞ | Undefined |

---

## Inverse Trigonometric Tests

### atan(x) and atan2(y, x)

| Function | Input(s) | Expected | Notes |
|----------|----------|----------|-------|
| `atan(0)` | 0 | `0` | |
| `atan(1)` | 1 | `π/4 = 0.785398...` | |
| `atan(-1)` | -1 | `-π/4` | |
| `atan(sqrt(3))` | √3 | `π/3` | |
| `atan(1/sqrt(3))` | 1/√3 | `π/6` | |
| `atan2(0, 1)` | (0, 1) | `0` | Positive x-axis |
| `atan2(1, 0)` | (1, 0) | `π/2` | Positive y-axis |
| `atan2(0, -1)` | (0, -1) | `π` | Negative x-axis |
| `atan2(-1, 0)` | (-1, 0) | `-π/2` | Negative y-axis |
| `atan2(1, 1)` | (1, 1) | `π/4` | First quadrant |
| `atan2(1, -1)` | (1, -1) | `3π/4` | Second quadrant |
| `atan2(-1, -1)` | (-1, -1) | `-3π/4` | Third quadrant |
| `atan2(-1, 1)` | (-1, 1) | `-π/4` | Fourth quadrant |
| `atan2(0, 0)` | (0, 0) | Error/NaN | Undefined |

### asin(x) and acos(x)

| x | asin(x) | acos(x) | Notes |
|---|---------|---------|-------|
| `0` | `0` | `π/2` | |
| `1` | `π/2` | `0` | |
| `-1` | `-π/2` | `π` | |
| `0.5` | `π/6` | `π/3` | |
| `sqrt(2)/2` | `π/4` | `π/4` | |
| `sqrt(3)/2` | `π/3` | `π/6` | |
| `1.1` | Error/NaN | Error/NaN | Out of domain |
| `-1.1` | Error/NaN | Error/NaN | Out of domain |

---

## Hyperbolic Tests

### sinh(x) and cosh(x)

| x | sinh(x) | cosh(x) | Notes |
|---|---------|---------|-------|
| `0` | `0` | `1` | Exact |
| `1` | `1.175201193643801...` | `1.543080634815243...` | |
| `-1` | `-1.175201193643801...` | `1.543080634815243...` | cosh is even |
| `log(2)` | `0.75` | `1.25` | sinh(log(2)) = (2-1/2)/2 |
| `10` | `11013.23287470339...` | `11013.23292010332...` | Large |
| `0.001` | `0.001000000166666...` | `1.0000005` | Small (Taylor) |

### tanh(x)

| x | tanh(x) | Notes |
|---|---------|-------|
| `0` | `0` | |
| `1` | `0.761594155955764...` | |
| `-1` | `-0.761594155955764...` | |
| `10` | `0.999999999586...` | Approaches 1 |
| `-10` | `-0.999999999586...` | Approaches -1 |
| `100` | `1.0 - tiny` | Very close to 1 |

---

## Inverse Hyperbolic Tests

### asinh(x), acosh(x), atanh(x)

| Function | Input | Expected | Notes |
|----------|-------|----------|-------|
| `asinh(0)` | 0 | `0` | |
| `asinh(1)` | 1 | `0.881373587019543...` | |
| `asinh(sinh(5))` | sinh(5) | `5` | Inverse test |
| `acosh(1)` | 1 | `0` | |
| `acosh(2)` | 2 | `1.316957896924816...` | |
| `acosh(cosh(3))` | cosh(3) | `3` | Inverse test |
| `acosh(0.5)` | 0.5 | Error/NaN | Domain error |
| `atanh(0)` | 0 | `0` | |
| `atanh(0.5)` | 0.5 | `0.549306144334054...` | |
| `atanh(tanh(2))` | tanh(2) | `2` | Inverse test |
| `atanh(1)` | 1 | Error/NaN | +∞ |
| `atanh(-1)` | -1 | Error/NaN | -∞ |
| `atanh(1.5)` | 1.5 | Error/NaN | Domain error |

---

## Identity Tests

These tests verify mathematical identities to full precision.

### Pythagorean Identities

| Identity | Test | Expected |
|----------|------|----------|
| `sin²(x) + cos²(x) = 1` | Random x ∈ [-10, 10] | `1` |
| `cosh²(x) - sinh²(x) = 1` | Random x ∈ [-10, 10] | `1` |
| `1 + tan²(x) = sec²(x)` | x where cos(x) ≠ 0 | Equal |
| `1 + cot²(x) = csc²(x)` | x where sin(x) ≠ 0 | Equal |

### Inverse Function Identities

| Identity | Domain | Notes |
|----------|--------|-------|
| `sin(asin(x)) = x` | \|x\| ≤ 1 | |
| `asin(sin(x)) = x` | \|x\| ≤ π/2 | Restricted domain |
| `cos(acos(x)) = x` | \|x\| ≤ 1 | |
| `exp(log(x)) = x` | x > 0 | |
| `log(exp(x)) = x` | All x | |
| `sinh(asinh(x)) = x` | All x | |
| `cosh(acosh(x)) = x` | x ≥ 1 | |
| `tanh(atanh(x)) = x` | \|x\| < 1 | |

### Addition Formulas

| Formula | Test |
|---------|------|
| `sin(a+b) = sin(a)cos(b) + cos(a)sin(b)` | Random a, b |
| `cos(a+b) = cos(a)cos(b) - sin(a)sin(b)` | Random a, b |
| `exp(a+b) = exp(a) * exp(b)` | Random a, b |
| `log(a*b) = log(a) + log(b)` | a, b > 0 |

### Double/Half Angle Formulas

| Formula | Test |
|---------|------|
| `sin(2x) = 2*sin(x)*cos(x)` | Random x |
| `cos(2x) = cos²(x) - sin²(x)` | Random x |
| `sin(x/2) = ±√((1-cos(x))/2)` | Random x |

---

## Edge Cases and Error Handling

### Domain Errors (should return NaN or error state)

| Function | Invalid Input | Expected |
|----------|---------------|----------|
| `log(0)` | 0 | NaN |
| `log(-5)` | -5 | NaN |
| `sqrt(-1)` | -1 | NaN (for real sqrt) |
| `asin(2)` | 2 | NaN |
| `acos(-2)` | -2 | NaN |
| `acosh(0.5)` | 0.5 | NaN |
| `atanh(1)` | 1 | NaN or ±∞ indication |
| `atan2(0, 0)` | (0, 0) | NaN |
| `pow(-2, 0.5)` | (-2, 0.5) | NaN |

### Overflow/Underflow

| Function | Input | Expected Behavior |
|----------|-------|-------------------|
| `exp(1e10)` | Very large | Overflow → error or max value |
| `exp(-1e10)` | Very negative | Underflow → 0 |
| `log(1e-1000000)` | Very small | Large negative value |
| `sin(1e100)` | Very large | Compute if possible, warn about precision |

### Zero and Sign Handling

| Test | Expected |
|------|----------|
| `sin(-0) = -0` | Preserve sign of zero if supported |
| `cos(-0) = 1` | |
| `atan2(-0, 1) = -0` | |
| `atan2(0, -1) = π` | |
| `atan2(-0, -1) = -π` | |

---

## Precision Validation

### Varying Precision Levels

Test each function at multiple precision levels:

| Precision | Words | Approx. Digits |
|-----------|-------|----------------|
| Low | 4 | ~70 |
| Medium | 10 | ~180 |
| High | 50 | ~900 |
| Very High | 200 | ~3600 |

### Convergence Tests

For iterative algorithms (Newton-Raphson, series):
- Verify iteration count is O(log(precision))
- Verify final precision matches requested precision

### Round-Trip Tests

```
x → f(x) → f⁻¹(f(x)) ≈ x
```

With tolerance appropriate for precision level.

---

## Cross-Validation Against Fortran MPFUN

### Generate Reference Values

Create a Fortran program to output reference values:

```fortran
program gen_reference
    use mpmodule
    implicit none
    type(mp_real) :: x, result
    character(len=1000) :: str
    
    ! Initialize to high precision
    call mpinit(200)  ! 200 words ≈ 3600 digits
    
    ! exp(1)
    x = mpreal('1.0')
    result = exp(x)
    call mpwrite(6, 100, 80, result)
    
    ! log(2)  
    x = mpreal('2.0')
    result = log(x)
    call mpwrite(6, 100, 80, result)
    
    ! sin(1)
    x = mpreal('1.0')
    result = sin(x)
    call mpwrite(6, 100, 80, result)
    
    ! ... more test cases
end program
```

### Comparison Framework

```cpp
// Pseudocode for Agent 3 test infrastructure
void validate_against_fortran(const char* function_name, 
                               const MPFloat& input,
                               const char* expected_str) {
    MPFloat expected = from_string(expected_str);
    MPFloat computed = call_function(function_name, input);
    
    // Allow small relative error
    double rel_error = abs(computed - expected) / abs(expected);
    ASSERT(rel_error < 1e-(precision_digits - 5));
}
```

---

## Test Categories Summary

| Category | Purpose | Count (approx) |
|----------|---------|----------------|
| Basic values | Known mathematical constants | 50+ |
| Special values | 0, 1, ±∞, NaN handling | 30+ |
| Domain errors | Invalid input handling | 20+ |
| Identity tests | Mathematical invariants | 30+ |
| Precision tests | Accuracy at various precisions | 20+ |
| Cross-validation | Compare with Fortran MPFUN | 50+ |
| Stress tests | Large arguments, edge cases | 20+ |

**Total estimated test cases:** 200+

---

## References

1. NIST Digital Library of Mathematical Functions (DLMF) - https://dlmf.nist.gov/
2. Abramowitz & Stegun - Handbook of Mathematical Functions
3. MPFUN2020 test suite (`testmpfun.f90`)
4. IEEE 754-2019 - Floating-point arithmetic standard (for edge case behavior)
