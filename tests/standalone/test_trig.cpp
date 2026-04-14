/**
 * @file test_trig.cpp
 * @brief Standalone test for MPFloat trigonometric functions
 * 
 * Compile with:
 *   clang++ -std=c++17 -I../../include -o test_trig test_trig.cpp
 * 
 * Or with more precision:
 *   clang++ -std=c++17 -O2 -I../../include -o test_trig test_trig.cpp
 */

#include <cstdio>
#include <cstdint>
#include <cmath>
#include <cstring>

// Include MPFloat headers
#include "mpfun/core/representation.hpp"
#include "mpfun/core/add.hpp"
#include "mpfun/core/mul.hpp"
#include "mpfun/core/div.hpp"
#include "mpfun/core/sqrt.hpp"
// Include mp_float.hpp before trig.hpp since trig uses MPFloat class
#include "mpfun/mp_float.hpp"
#include "mpfun/transcendental/trig.hpp"

using namespace mpfun;
using namespace mpfun::detail;

constexpr int WORDS = 16;  // ~280 decimal digits
using MP = MPFloat<WORDS>;

// Reference π value for testing (more digits than double can hold)
constexpr double PI = 3.14159265358979323846;
constexpr double PI_2 = PI / 2.0;
constexpr double PI_4 = PI / 4.0;
constexpr double PI_6 = PI / 6.0;
constexpr double PI_3 = PI / 3.0;

// Helper to print an MPFloat value with name
void print_mp(const char* name, const MP& x) {
    double d = x.to_double();
    int s = x.sign();
    int len = x.mantissa_length();
    int64_t exp = x.exponent();
    
    printf("  %s: %.15g (sign=%d, len=%d, exp=%lld)\n", 
           name, d, s, len, (long long)exp);
}

// Test result tracking
int tests_passed = 0;
int tests_failed = 0;

bool check_close(const char* test_name, const MP& result, double expected, double rel_tol = 1e-10) {
    double actual = result.to_double();
    double abs_diff = std::fabs(actual - expected);
    double rel_diff = (expected != 0.0) ? abs_diff / std::fabs(expected) : abs_diff;
    
    bool passed = (rel_diff < rel_tol) || (abs_diff < 1e-15);
    
    if (passed) {
        printf("✓ %s: %.15g ≈ %.15g (rel_err=%.2e)\n", 
               test_name, actual, expected, rel_diff);
        tests_passed++;
    } else {
        printf("✗ %s: %.15g != %.15g (rel_err=%.2e, tol=%.2e)\n", 
               test_name, actual, expected, rel_diff, rel_tol);
        tests_failed++;
    }
    return passed;
}

bool check_identity(const char* test_name, const MP& result, double expected, double abs_tol = 1e-12) {
    double actual = result.to_double();
    double abs_diff = std::fabs(actual - expected);
    
    bool passed = (abs_diff < abs_tol);
    
    if (passed) {
        printf("✓ %s: %.15g ≈ %.15g (abs_err=%.2e)\n", 
               test_name, actual, expected, abs_diff);
        tests_passed++;
    } else {
        printf("✗ %s: %.15g != %.15g (abs_err=%.2e, tol=%.2e)\n", 
               test_name, actual, expected, abs_diff, abs_tol);
        tests_failed++;
    }
    return passed;
}

// =============================================================================
// Test: sin and cos basic values
// =============================================================================
void test_sincos_basic() {
    printf("\n=== Test: sin/cos basic values ===\n");
    
    // sin(0) = 0
    {
        MP x(0.0);
        MP s = sin(x);
        check_identity("sin(0) = 0", s, 0.0);
    }
    
    // cos(0) = 1
    {
        MP x(0.0);
        MP c = cos(x);
        check_identity("cos(0) = 1", c, 1.0);
    }
    
    // sin(π/6) = 0.5
    {
        MP x(PI_6);
        MP s = sin(x);
        check_close("sin(π/6) = 0.5", s, 0.5);
    }
    
    // cos(π/3) = 0.5
    {
        MP x(PI_3);
        MP c = cos(x);
        check_close("cos(π/3) = 0.5", c, 0.5);
    }
    
    // sin(π/4) = cos(π/4) = √2/2
    {
        MP x(PI_4);
        MP s = sin(x);
        MP c = cos(x);
        double sqrt2_2 = std::sqrt(2.0) / 2.0;
        check_close("sin(π/4) = √2/2", s, sqrt2_2);
        check_close("cos(π/4) = √2/2", c, sqrt2_2);
    }
    
    // sin(π/2) = 1
    {
        MP x(PI_2);
        MP s = sin(x);
        check_close("sin(π/2) = 1", s, 1.0);
    }
    
    // cos(π/2) = 0
    {
        MP x(PI_2);
        MP c = cos(x);
        check_identity("cos(π/2) ≈ 0", c, 0.0, 1e-10);
    }
    
    // sin(π) ≈ 0
    {
        MP x(PI);
        MP s = sin(x);
        check_identity("sin(π) ≈ 0", s, 0.0, 1e-10);
    }
    
    // cos(π) = -1
    {
        MP x(PI);
        MP c = cos(x);
        check_close("cos(π) = -1", c, -1.0);
    }
}

// =============================================================================
// Test: sin²(x) + cos²(x) = 1 identity
// =============================================================================
void test_pythagorean_identity() {
    printf("\n=== Test: sin²(x) + cos²(x) = 1 ===\n");
    
    double test_angles[] = {0.0, 0.1, 0.5, 1.0, PI_4, PI_3, PI_2, PI, 2*PI, 3.0, 10.0, -1.0, -PI};
    
    for (double angle : test_angles) {
        MP x(angle);
        MP s, c;
        sincos(x, s, c);
        
        MP s_sq = s * s;
        MP c_sq = c * c;
        MP sum = s_sq + c_sq;
        
        char buf[64];
        snprintf(buf, sizeof(buf), "sin²(%.2f) + cos²(%.2f) = 1", angle, angle);
        check_identity(buf, sum, 1.0, 1e-10);
    }
}

// =============================================================================
// Test: tan(x) = sin(x) / cos(x)
// =============================================================================
void test_tan() {
    printf("\n=== Test: tan(x) ===\n");
    
    // tan(0) = 0
    {
        MP x(0.0);
        MP t = tan(x);
        check_identity("tan(0) = 0", t, 0.0);
    }
    
    // tan(π/4) = 1
    {
        MP x(PI_4);
        MP t = tan(x);
        check_close("tan(π/4) = 1", t, 1.0);
    }
    
    // tan(π/6) = 1/√3
    {
        MP x(PI_6);
        MP t = tan(x);
        check_close("tan(π/6) = 1/√3", t, 1.0 / std::sqrt(3.0));
    }
    
    // tan(π/3) = √3
    {
        MP x(PI_3);
        MP t = tan(x);
        check_close("tan(π/3) = √3", t, std::sqrt(3.0));
    }
    
    // tan(-π/4) = -1
    {
        MP x(-PI_4);
        MP t = tan(x);
        check_close("tan(-π/4) = -1", t, -1.0);
    }
}

// =============================================================================
// Test: atan(x)
// =============================================================================
void test_atan() {
    printf("\n=== Test: atan(x) ===\n");
    
    // atan(0) = 0
    {
        MP x(0.0);
        MP a = atan(x);
        check_identity("atan(0) = 0", a, 0.0);
    }
    
    // atan(1) = π/4
    {
        MP x(1.0);
        MP a = atan(x);
        check_close("atan(1) = π/4", a, PI_4);
    }
    
    // atan(-1) = -π/4
    {
        MP x(-1.0);
        MP a = atan(x);
        check_close("atan(-1) = -π/4", a, -PI_4);
    }
    
    // atan(√3) = π/3
    {
        MP x(std::sqrt(3.0));
        MP a = atan(x);
        check_close("atan(√3) = π/3", a, PI_3);
    }
    
    // atan(1/√3) = π/6
    {
        MP x(1.0 / std::sqrt(3.0));
        MP a = atan(x);
        check_close("atan(1/√3) = π/6", a, PI_6);
    }
    
    // atan(large) ≈ π/2
    {
        MP x(1000.0);
        MP a = atan(x);
        check_close("atan(1000) ≈ π/2", a, PI_2, 1e-3);
    }
}

// =============================================================================
// Test: atan2(y, x)
// =============================================================================
void test_atan2() {
    printf("\n=== Test: atan2(y, x) ===\n");
    
    // atan2(0, 1) = 0
    {
        MP y(0.0), x(1.0);
        MP a = atan2(y, x);
        check_identity("atan2(0, 1) = 0", a, 0.0);
    }
    
    // atan2(1, 0) = π/2
    {
        MP y(1.0), x(0.0);
        MP a = atan2(y, x);
        check_close("atan2(1, 0) = π/2", a, PI_2);
    }
    
    // atan2(-1, 0) = -π/2
    {
        MP y(-1.0), x(0.0);
        MP a = atan2(y, x);
        check_close("atan2(-1, 0) = -π/2", a, -PI_2);
    }
    
    // atan2(0, -1) = π
    {
        MP y(0.0), x(-1.0);
        MP a = atan2(y, x);
        check_close("atan2(0, -1) = π", a, PI);
    }
    
    // atan2(1, 1) = π/4
    {
        MP y(1.0), x(1.0);
        MP a = atan2(y, x);
        check_close("atan2(1, 1) = π/4", a, PI_4);
    }
    
    // atan2(-1, -1) = -3π/4
    {
        MP y(-1.0), x(-1.0);
        MP a = atan2(y, x);
        check_close("atan2(-1, -1) = -3π/4", a, -3.0 * PI_4);
    }
}

// =============================================================================
// Test: asin(x)
// =============================================================================
void test_asin() {
    printf("\n=== Test: asin(x) ===\n");
    
    // asin(0) = 0
    {
        MP x(0.0);
        MP a = asin(x);
        check_identity("asin(0) = 0", a, 0.0);
    }
    
    // asin(1) = π/2
    {
        MP x(1.0);
        MP a = asin(x);
        check_close("asin(1) = π/2", a, PI_2);
    }
    
    // asin(-1) = -π/2
    {
        MP x(-1.0);
        MP a = asin(x);
        check_close("asin(-1) = -π/2", a, -PI_2);
    }
    
    // asin(0.5) = π/6
    {
        MP x(0.5);
        MP a = asin(x);
        check_close("asin(0.5) = π/6", a, PI_6);
    }
    
    // asin(√2/2) = π/4
    {
        MP x(std::sqrt(2.0) / 2.0);
        MP a = asin(x);
        check_close("asin(√2/2) = π/4", a, PI_4);
    }
    
    // asin(√3/2) = π/3
    {
        MP x(std::sqrt(3.0) / 2.0);
        MP a = asin(x);
        check_close("asin(√3/2) = π/3", a, PI_3);
    }
}

// =============================================================================
// Test: acos(x)
// =============================================================================
void test_acos() {
    printf("\n=== Test: acos(x) ===\n");
    
    // acos(1) = 0
    {
        MP x(1.0);
        MP a = acos(x);
        check_identity("acos(1) = 0", a, 0.0, 1e-10);
    }
    
    // acos(0) = π/2
    {
        MP x(0.0);
        MP a = acos(x);
        check_close("acos(0) = π/2", a, PI_2);
    }
    
    // acos(-1) = π
    {
        MP x(-1.0);
        MP a = acos(x);
        check_close("acos(-1) = π", a, PI);
    }
    
    // acos(0.5) = π/3
    {
        MP x(0.5);
        MP a = acos(x);
        check_close("acos(0.5) = π/3", a, PI_3);
    }
    
    // acos(√2/2) = π/4
    {
        MP x(std::sqrt(2.0) / 2.0);
        MP a = acos(x);
        check_close("acos(√2/2) = π/4", a, PI_4);
    }
    
    // acos(√3/2) = π/6
    {
        MP x(std::sqrt(3.0) / 2.0);
        MP a = acos(x);
        check_close("acos(√3/2) = π/6", a, PI_6);
    }
}

// =============================================================================
// Test: Inverse function roundtrips
// =============================================================================
void test_inverse_roundtrip() {
    printf("\n=== Test: Inverse function roundtrips ===\n");
    
    double test_values[] = {0.0, 0.1, 0.25, 0.5, 0.75, 0.9, 0.99};
    
    for (double v : test_values) {
        // sin(asin(x)) = x
        {
            MP x(v);
            MP a = asin(x);
            MP result = sin(a);
            char buf[64];
            snprintf(buf, sizeof(buf), "sin(asin(%.2f)) = %.2f", v, v);
            check_close(buf, result, v);
        }
        
        // cos(acos(x)) = x
        {
            MP x(v);
            MP a = acos(x);
            MP result = cos(a);
            char buf[64];
            snprintf(buf, sizeof(buf), "cos(acos(%.2f)) = %.2f", v, v);
            check_close(buf, result, v);
        }
    }
    
    // tan(atan(x)) = x for various values
    double tan_values[] = {0.0, 0.1, 0.5, 1.0, 2.0, 5.0, 10.0, -1.0, -0.5};
    for (double v : tan_values) {
        MP x(v);
        MP a = atan(x);
        MP result = tan(a);
        char buf[64];
        snprintf(buf, sizeof(buf), "tan(atan(%.2f)) = %.2f", v, v);
        check_close(buf, result, v);
    }
}

// =============================================================================
// Test: Large argument reduction
// =============================================================================
void test_large_arguments() {
    printf("\n=== Test: Large argument reduction ===\n");
    
    // sin(2πn) = 0 for integer n
    for (int n = 1; n <= 5; ++n) {
        MP x(2.0 * PI * n);
        MP s = sin(x);
        char buf[64];
        snprintf(buf, sizeof(buf), "sin(2π×%d) ≈ 0", n);
        check_identity(buf, s, 0.0, 1e-8);
    }
    
    // cos(2πn) = 1 for integer n
    for (int n = 1; n <= 5; ++n) {
        MP x(2.0 * PI * n);
        MP c = cos(x);
        char buf[64];
        snprintf(buf, sizeof(buf), "cos(2π×%d) = 1", n);
        check_close(buf, c, 1.0, 1e-8);
    }
    
    // sin(100) - compare with standard library
    {
        MP x(100.0);
        MP s = sin(x);
        check_close("sin(100) vs std::sin", s, std::sin(100.0), 1e-8);
    }
    
    // cos(100) - compare with standard library
    {
        MP x(100.0);
        MP c = cos(x);
        check_close("cos(100) vs std::cos", c, std::cos(100.0), 1e-8);
    }
}

// =============================================================================
// Test: sincos efficiency check
// =============================================================================
void test_sincos_direct() {
    printf("\n=== Test: sincos() direct call ===\n");
    
    MP x(1.0);
    MP s, c;
    
    sincos(x, s, c);
    
    check_close("sincos: sin(1) vs std::sin", s, std::sin(1.0));
    check_close("sincos: cos(1) vs std::cos", c, std::cos(1.0));
    
    // Verify Pythagorean identity
    MP sum = s * s + c * c;
    check_identity("sincos: sin²(1) + cos²(1) = 1", sum, 1.0, 1e-10);
}

// =============================================================================
// Main
// =============================================================================
int main() {
    printf("MPFloat Trigonometric Functions Test\n");
    printf("====================================\n");
    printf("Using %d words = ~%d decimal digits\n\n", WORDS, WORDS * 18);
    
    test_sincos_basic();
    test_pythagorean_identity();
    test_tan();
    test_atan();
    test_atan2();
    test_asin();
    test_acos();
    test_inverse_roundtrip();
    test_large_arguments();
    test_sincos_direct();
    
    printf("\n====================================\n");
    printf("Results: %d passed, %d failed\n", tests_passed, tests_failed);
    
    return (tests_failed == 0) ? 0 : 1;
}
