/**
 * @file test_hyperbolic.cpp
 * @brief Standalone tests for hyperbolic functions
 * 
 * Tests the hyperbolic functions implementation:
 * - sinh, cosh, tanh
 * - asinh, acosh, atanh
 * - sinhcosh (joint computation)
 * 
 * Build: clang++ -std=c++17 -I../../include -o test_hyperbolic test_hyperbolic.cpp
 * Run: ./test_hyperbolic
 */

#include <iostream>
#include <iomanip>
#include <cmath>
#include <string>

// Include MPFloat and hyperbolic functions
#include "mpfun/mp_float.hpp"
#include "mpfun/core/add.hpp"
#include "mpfun/core/mul.hpp"
#include "mpfun/core/div.hpp"
#include "mpfun/core/sqrt.hpp"
#include "mpfun/core/compare.hpp"
#include "mpfun/transcendental/hyperbolic.hpp"

using namespace mpfun;

// Test precision: ~100 decimal digits
constexpr int WORDS = 6;
using MP = MPFloat<WORDS>;

// Helper to check if two values are approximately equal
bool approx_equal(const MP& a, const MP& b, double rel_tol = 1e-15) {
    MP diff = a - b;
    MP abs_diff = abs(diff);
    MP max_val = abs(a);
    MP abs_b = abs(b);
    if (abs_b > max_val) max_val = abs_b;
    
    // Handle near-zero case
    if (max_val.to_double() < 1e-100) {
        return abs_diff.to_double() < rel_tol;
    }
    
    return (abs_diff / max_val).to_double() < rel_tol;
}

bool approx_equal_double(double a, double b, double rel_tol = 1e-12) {
    if (std::abs(a) < 1e-100 && std::abs(b) < 1e-100) {
        return std::abs(a - b) < rel_tol;
    }
    double max_val = std::max(std::abs(a), std::abs(b));
    return std::abs(a - b) / max_val < rel_tol;
}

int tests_passed = 0;
int tests_failed = 0;

void test(bool condition, const std::string& name) {
    if (condition) {
        std::cout << "  ✓ " << name << std::endl;
        tests_passed++;
    } else {
        std::cout << "  ✗ " << name << " FAILED" << std::endl;
        tests_failed++;
    }
}

void test_section(const std::string& name) {
    std::cout << "\n" << name << ":" << std::endl;
}

// =============================================================================
// Test Basic Hyperbolic Functions
// =============================================================================

void test_sinh_basic() {
    test_section("sinh() basic tests");
    
    // sinh(0) = 0
    MP zero(0);
    MP sinh_zero = sinh(zero);
    test(std::abs(sinh_zero.to_double()) < 1e-50, "sinh(0) = 0");
    
    // sinh is odd: sinh(-x) = -sinh(x)
    MP x(1.5);
    MP sinh_x = sinh(x);
    MP neg_x = -x;
    MP sinh_neg_x = sinh(neg_x);
    test(approx_equal(sinh_neg_x, -sinh_x), "sinh(-x) = -sinh(x)");
    
    // Compare with std::sinh for moderate values
    double val = 0.5;
    MP mp_val(val);
    double mp_sinh = sinh(mp_val).to_double();
    double std_sinh = std::sinh(val);
    test(approx_equal_double(mp_sinh, std_sinh), "sinh(0.5) matches std::sinh");
    
    val = 2.0;
    mp_val = MP(val);
    mp_sinh = sinh(mp_val).to_double();
    std_sinh = std::sinh(val);
    test(approx_equal_double(mp_sinh, std_sinh), "sinh(2.0) matches std::sinh");
}

void test_cosh_basic() {
    test_section("cosh() basic tests");
    
    // cosh(0) = 1
    MP zero(0);
    MP cosh_zero = cosh(zero);
    test(approx_equal_double(cosh_zero.to_double(), 1.0), "cosh(0) = 1");
    
    // cosh is even: cosh(-x) = cosh(x)
    MP x(1.5);
    MP cosh_x = cosh(x);
    MP neg_x = -x;
    MP cosh_neg_x = cosh(neg_x);
    test(approx_equal(cosh_neg_x, cosh_x), "cosh(-x) = cosh(x)");
    
    // cosh(x) >= 1 for all x
    MP cosh_5 = cosh(MP(5.0));
    test(cosh_5.to_double() >= 1.0, "cosh(5) >= 1");
    
    // Compare with std::cosh
    double val = 1.0;
    MP mp_val(val);
    double mp_cosh = cosh(mp_val).to_double();
    double std_cosh = std::cosh(val);
    test(approx_equal_double(mp_cosh, std_cosh), "cosh(1.0) matches std::cosh");
}

void test_fundamental_identity() {
    test_section("Fundamental identity: cosh²(x) - sinh²(x) = 1");
    
    // Test for several values
    double test_vals[] = {0.0, 0.001, 0.1, 0.5, 1.0, 2.0, 3.0, 5.0};
    
    for (double val : test_vals) {
        MP x(val);
        MP sinh_x, cosh_x;
        sinhcosh(x, sinh_x, cosh_x);
        
        MP cosh_sq = cosh_x * cosh_x;
        MP sinh_sq = sinh_x * sinh_x;
        MP diff = cosh_sq - sinh_sq;
        
        std::string name = "cosh²(" + std::to_string(val) + ") - sinh²(" + 
                           std::to_string(val) + ") = 1";
        test(approx_equal_double(diff.to_double(), 1.0, 1e-10), name);
    }
    
    // Test with negative values
    MP x(-2.5);
    MP sinh_x, cosh_x;
    sinhcosh(x, sinh_x, cosh_x);
    MP diff = cosh_x * cosh_x - sinh_x * sinh_x;
    test(approx_equal_double(diff.to_double(), 1.0, 1e-10), 
         "cosh²(-2.5) - sinh²(-2.5) = 1");
}

void test_tanh_basic() {
    test_section("tanh() basic tests");
    
    // tanh(0) = 0
    MP zero(0);
    MP tanh_zero = tanh(zero);
    test(std::abs(tanh_zero.to_double()) < 1e-50, "tanh(0) = 0");
    
    // tanh is odd: tanh(-x) = -tanh(x)
    MP x(1.0);
    MP tanh_x = tanh(x);
    MP tanh_neg_x = tanh(-x);
    test(approx_equal(tanh_neg_x, -tanh_x), "tanh(-x) = -tanh(x)");
    
    // |tanh(x)| < 1 for all finite x
    double large_val = 5.0;
    MP tanh_large = tanh(MP(large_val));
    test(std::abs(tanh_large.to_double()) < 1.0, "|tanh(5)| < 1");
    
    // tanh approaches ±1 for large |x|
    test(tanh_large.to_double() > 0.9999, "tanh(5) > 0.9999");
    
    // Compare with std::tanh
    double val = 0.5;
    double mp_tanh = tanh(MP(val)).to_double();
    double std_tanh = std::tanh(val);
    test(approx_equal_double(mp_tanh, std_tanh), "tanh(0.5) matches std::tanh");
}

void test_sinhcosh_joint() {
    test_section("sinhcosh() joint computation");
    
    MP x(1.5);
    MP sinh_joint, cosh_joint;
    sinhcosh(x, sinh_joint, cosh_joint);
    
    // Compare with separate computations
    MP sinh_sep = sinh(x);
    MP cosh_sep = cosh(x);
    
    test(approx_equal(sinh_joint, sinh_sep), "sinhcosh sinh == separate sinh");
    test(approx_equal(cosh_joint, cosh_sep), "sinhcosh cosh == separate cosh");
}

// =============================================================================
// Test Inverse Hyperbolic Functions
// =============================================================================

void test_asinh_basic() {
    test_section("asinh() basic tests");
    
    // asinh(0) = 0
    MP zero(0);
    MP asinh_zero = asinh(zero);
    test(std::abs(asinh_zero.to_double()) < 1e-50, "asinh(0) = 0");
    
    // asinh is odd: asinh(-x) = -asinh(x)
    MP x(2.0);
    MP asinh_x = asinh(x);
    MP asinh_neg_x = asinh(-x);
    test(approx_equal(asinh_neg_x, -asinh_x), "asinh(-x) = -asinh(x)");
    
    // Compare with std::asinh
    double val = 1.0;
    double mp_asinh = asinh(MP(val)).to_double();
    double std_asinh = std::asinh(val);
    test(approx_equal_double(mp_asinh, std_asinh), "asinh(1.0) matches std::asinh");
}

void test_asinh_inverse() {
    test_section("asinh(sinh(x)) = x identity");
    
    double test_vals[] = {0.0, 0.1, 0.5, 1.0, 2.0, 3.0};
    
    for (double val : test_vals) {
        MP x(val);
        MP sinh_x = sinh(x);
        MP asinh_sinh_x = asinh(sinh_x);
        
        std::string name = "asinh(sinh(" + std::to_string(val) + ")) = " + 
                           std::to_string(val);
        test(approx_equal_double(asinh_sinh_x.to_double(), val, 1e-10), name);
    }
    
    // Negative value
    MP x(-1.5);
    MP result = asinh(sinh(x));
    test(approx_equal_double(result.to_double(), -1.5, 1e-10), 
         "asinh(sinh(-1.5)) = -1.5");
}

void test_acosh_basic() {
    test_section("acosh() basic tests");
    
    // acosh(1) = 0
    MP one(1);
    MP acosh_one = acosh(one);
    test(std::abs(acosh_one.to_double()) < 1e-50, "acosh(1) = 0");
    
    // Compare with std::acosh
    double val = 2.0;
    double mp_acosh = acosh(MP(val)).to_double();
    double std_acosh = std::acosh(val);
    test(approx_equal_double(mp_acosh, std_acosh), "acosh(2.0) matches std::acosh");
    
    val = 1.5;
    mp_acosh = acosh(MP(val)).to_double();
    std_acosh = std::acosh(val);
    test(approx_equal_double(mp_acosh, std_acosh), "acosh(1.5) matches std::acosh");
}

void test_acosh_inverse() {
    test_section("acosh(cosh(x)) = x for x > 0");
    
    // acosh(cosh(x)) = |x| since cosh is even and acosh returns non-negative branch
    double test_vals[] = {0.1, 0.5, 1.0, 2.0, 3.0};
    
    for (double val : test_vals) {
        MP x(val);
        MP cosh_x = cosh(x);
        MP acosh_cosh_x = acosh(cosh_x);
        
        std::string name = "acosh(cosh(" + std::to_string(val) + ")) = " + 
                           std::to_string(val);
        test(approx_equal_double(acosh_cosh_x.to_double(), val, 1e-10), name);
    }
}

void test_acosh_domain() {
    test_section("acosh() domain check");
    
    // acosh(x) for x < 1 should return NaN
    MP less_than_one(0.5);
    MP result = acosh(less_than_one);
    test(result.is_nan(), "acosh(0.5) returns NaN");
}

void test_atanh_basic() {
    test_section("atanh() basic tests");
    
    // atanh(0) = 0
    MP zero(0);
    MP atanh_zero = atanh(zero);
    test(std::abs(atanh_zero.to_double()) < 1e-50, "atanh(0) = 0");
    
    // atanh is odd: atanh(-x) = -atanh(x)
    MP x(0.5);
    MP atanh_x = atanh(x);
    MP atanh_neg_x = atanh(-x);
    test(approx_equal(atanh_neg_x, -atanh_x), "atanh(-x) = -atanh(x)");
    
    // Compare with std::atanh
    double val = 0.5;
    double mp_atanh = atanh(MP(val)).to_double();
    double std_atanh = std::atanh(val);
    test(approx_equal_double(mp_atanh, std_atanh), "atanh(0.5) matches std::atanh");
    
    val = 0.9;
    mp_atanh = atanh(MP(val)).to_double();
    std_atanh = std::atanh(val);
    test(approx_equal_double(mp_atanh, std_atanh), "atanh(0.9) matches std::atanh");
}

void test_atanh_domain() {
    test_section("atanh() domain check");
    
    // atanh(x) for |x| >= 1 should return NaN
    MP one(1.0);
    MP result = atanh(one);
    test(result.is_nan(), "atanh(1.0) returns NaN");
    
    MP greater_than_one(1.5);
    result = atanh(greater_than_one);
    test(result.is_nan(), "atanh(1.5) returns NaN");
}

// =============================================================================
// Test Small Argument Behavior (Taylor Series Path)
// =============================================================================

void test_small_arguments() {
    test_section("Small argument handling (Taylor series)");
    
    // For very small x, sinh(x) ≈ x
    MP tiny(1e-10);
    MP sinh_tiny = sinh(tiny);
    test(approx_equal_double(sinh_tiny.to_double(), 1e-10, 1e-5), 
         "sinh(1e-10) ≈ 1e-10");
    
    // For very small x, cosh(x) ≈ 1
    MP cosh_tiny = cosh(tiny);
    test(approx_equal_double(cosh_tiny.to_double(), 1.0, 1e-15), 
         "cosh(1e-10) ≈ 1");
    
    // For very small x, tanh(x) ≈ x
    MP tanh_tiny = tanh(tiny);
    test(approx_equal_double(tanh_tiny.to_double(), 1e-10, 1e-5), 
         "tanh(1e-10) ≈ 1e-10");
    
    // Identity should still hold
    MP sinh_val, cosh_val;
    sinhcosh(tiny, sinh_val, cosh_val);
    MP diff = cosh_val * cosh_val - sinh_val * sinh_val;
    test(approx_equal_double(diff.to_double(), 1.0, 1e-15), 
         "cosh²(1e-10) - sinh²(1e-10) = 1");
}

// =============================================================================
// Test Large Argument Behavior (Exponential Path)
// =============================================================================

void test_large_arguments() {
    test_section("Large argument handling (exponential)");
    
    // For large x, sinh(x) ≈ cosh(x) ≈ exp(x)/2
    MP large(10.0);
    MP sinh_large = sinh(large);
    MP cosh_large = cosh(large);
    
    double exp_10_half = std::exp(10.0) / 2.0;
    test(approx_equal_double(sinh_large.to_double(), exp_10_half, 1e-5), 
         "sinh(10) ≈ exp(10)/2");
    test(approx_equal_double(cosh_large.to_double(), exp_10_half, 1e-5), 
         "cosh(10) ≈ exp(10)/2");
    
    // tanh(large) ≈ 1
    MP tanh_large = tanh(large);
    test(tanh_large.to_double() > 0.9999999, "tanh(10) ≈ 1");
    
    // Identity should still hold
    MP diff = cosh_large * cosh_large - sinh_large * sinh_large;
    test(approx_equal_double(diff.to_double(), 1.0, 1e-8), 
         "cosh²(10) - sinh²(10) = 1");
}

// =============================================================================
// Main
// =============================================================================

int main() {
    std::cout << "=== MPFloat Hyperbolic Functions Test Suite ===" << std::endl;
    std::cout << "Precision: " << WORDS << " words (~" 
              << MP::approx_digits() << " decimal digits)" << std::endl;
    
    // Basic function tests
    test_sinh_basic();
    test_cosh_basic();
    test_fundamental_identity();
    test_tanh_basic();
    test_sinhcosh_joint();
    
    // Inverse function tests
    test_asinh_basic();
    test_asinh_inverse();
    test_acosh_basic();
    test_acosh_inverse();
    test_acosh_domain();
    test_atanh_basic();
    test_atanh_domain();
    
    // Edge case tests
    test_small_arguments();
    test_large_arguments();
    
    // Summary
    std::cout << "\n=== Summary ===" << std::endl;
    std::cout << "Tests passed: " << tests_passed << std::endl;
    std::cout << "Tests failed: " << tests_failed << std::endl;
    
    return tests_failed > 0 ? 1 : 0;
}
