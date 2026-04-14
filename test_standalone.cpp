/**
 * @file test_standalone.cpp
 * @brief Standalone compilation and basic sanity test
 * 
 * Compile with: g++ -std=c++17 -I include test_standalone.cpp -o test_standalone
 */

#include <iostream>
#include <iomanip>
#include <cmath>
#include "include/mpfun/mpfun.hpp"

using namespace mpfun;

// Helper to print test results
#define TEST(name, cond) do { \
    if (cond) { \
        std::cout << "[PASS] " << name << std::endl; \
        passed++; \
    } else { \
        std::cout << "[FAIL] " << name << std::endl; \
        failed++; \
    } \
} while(0)

int main() {
    int passed = 0, failed = 0;
    
    std::cout << "=== MPFloat Standalone Test ===" << std::endl;
    std::cout << std::setprecision(15);
    
    // Basic construction
    {
        MPFloat100 x;
        TEST("Default constructor creates zero", x.is_zero());
    }
    
    {
        MPFloat100 a(3.14159265358979);
        TEST("Double constructor", std::abs(a.to_double() - 3.14159265358979) < 1e-14);
    }
    
    // Addition
    {
        MPFloat100 a(100.0);
        MPFloat100 b(23.5);
        MPFloat100 c = a + b;
        TEST("Addition 100 + 23.5 = 123.5", std::abs(c.to_double() - 123.5) < 1e-10);
    }
    
    // Subtraction
    {
        MPFloat100 a(100.0);
        MPFloat100 b(23.5);
        MPFloat100 c = a - b;
        TEST("Subtraction 100 - 23.5 = 76.5", std::abs(c.to_double() - 76.5) < 1e-10);
    }
    
    // Multiplication
    {
        MPFloat100 a(12.0);
        MPFloat100 b(3.5);
        MPFloat100 c = a * b;
        TEST("Multiplication 12 * 3.5 = 42", std::abs(c.to_double() - 42.0) < 1e-10);
    }
    
    // Division
    {
        MPFloat100 a(100.0);
        MPFloat100 b(4.0);
        MPFloat100 c = a / b;
        TEST("Division 100 / 4 = 25", std::abs(c.to_double() - 25.0) < 1e-10);
    }
    
    {
        MPFloat100 one(1.0);
        MPFloat100 two(2.0);
        MPFloat100 half = one / two;
        TEST("Division 1 / 2 = 0.5", std::abs(half.to_double() - 0.5) < 1e-14);
    }
    
    {
        MPFloat100 one(1.0);
        MPFloat100 three(3.0);
        MPFloat100 third = one / three;
        TEST("Division 1 / 3", std::abs(third.to_double() - 1.0/3.0) < 1e-14);
        
        // Verify roundtrip
        MPFloat100 back = third * three;
        TEST("Roundtrip 3 * (1/3) ≈ 1", std::abs(back.to_double() - 1.0) < 1e-12);
    }
    
    // Division edge cases
    {
        MPFloat100 zero;
        MPFloat100 nonzero(42.0);
        MPFloat100 result = zero / nonzero;
        TEST("0 / x = 0", result.is_zero());
    }
    
    {
        MPFloat100 a(42.0);
        MPFloat100 zero;
        MPFloat100 result = a / zero;
        TEST("x / 0 = NaN", result.is_nan());
    }
    
    // Newton-Raphson convergence test for division
    {
        MPFloat100 a(7.0);
        MPFloat100 one(1.0);
        MPFloat100 inv_a = one / a;
        MPFloat100 product = a * inv_a;
        MPFloat100 error = product - one;
        double err_double = error.to_double();
        TEST("Newton-Raphson div convergence: 7 * (1/7) - 1", std::abs(err_double) < 1e-14);
    }
    
    // Square root
    {
        MPFloat100 four(4.0);
        MPFloat100 result = sqrt(four);
        TEST("sqrt(4) = 2", std::abs(result.to_double() - 2.0) < 1e-14);
    }
    
    {
        MPFloat100 nine(9.0);
        MPFloat100 result = sqrt(nine);
        TEST("sqrt(9) = 3", std::abs(result.to_double() - 3.0) < 1e-14);
    }
    
    {
        MPFloat100 two(2.0);
        MPFloat100 result = sqrt(two);
        TEST("sqrt(2)", std::abs(result.to_double() - std::sqrt(2.0)) < 1e-14);
        
        // Verify: result^2 ≈ 2
        MPFloat100 squared = result * result;
        TEST("sqrt(2)^2 ≈ 2", std::abs(squared.to_double() - 2.0) < 1e-12);
    }
    
    {
        MPFloat100 zero;
        MPFloat100 result = sqrt(zero);
        TEST("sqrt(0) = 0", result.is_zero());
    }
    
    {
        MPFloat100 neg(-4.0);
        MPFloat100 result = sqrt(neg);
        TEST("sqrt(-4) = NaN", result.is_nan());
    }
    
    // Combined test: Pythagorean theorem
    {
        MPFloat100 three(3.0);
        MPFloat100 four(4.0);
        MPFloat100 sum_squares = three * three + four * four;
        MPFloat100 hyp = sqrt(sum_squares);
        TEST("3^2 + 4^2 = 25, sqrt = 5", std::abs(hyp.to_double() - 5.0) < 1e-12);
    }
    
    // Division with exponent handling
    {
        MPFloat100 large(1e20);
        MPFloat100 small(1e10);
        MPFloat100 ratio = large / small;
        TEST("1e20 / 1e10 = 1e10", std::abs(ratio.to_double() / 1e10 - 1.0) < 1e-10);
    }
    
    std::cout << "\n=== Summary ===" << std::endl;
    std::cout << "Passed: " << passed << std::endl;
    std::cout << "Failed: " << failed << std::endl;
    
    return (failed == 0) ? 0 : 1;
}
