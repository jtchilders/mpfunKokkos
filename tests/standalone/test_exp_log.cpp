/**
 * @file test_exp_log.cpp
 * @brief Standalone tests for exp and log functions
 * 
 * Tests:
 * - exp(0) = 1
 * - exp(1) ≈ e
 * - log(1) = 0
 * - log(e) = 1
 * - exp(log(x)) = x
 * - log(exp(x)) = x
 * 
 * Compile with:
 *   g++ -std=c++17 -I../../include -O2 -o test_exp_log test_exp_log.cpp
 * 
 * Or with Kokkos:
 *   g++ -std=c++17 -I../../include -I$KOKKOS_PATH/include -O2 -o test_exp_log test_exp_log.cpp
 */

#include <iostream>
#include <iomanip>
#include <cmath>
#include <string>

// Include the full library
#include "mpfun/mp_float.hpp"
#include "mpfun/transcendental/exp_log.hpp"

using namespace mpfun;

// Reference values computed with high precision
constexpr double E_REF = 2.718281828459045235360287471352662497757;
constexpr double LN2_REF = 0.693147180559945309417232121458176568;
constexpr double LN10_REF = 2.302585092994045684017991454684364208;

// Test result tracking
int tests_passed = 0;
int tests_failed = 0;

void check(bool condition, const std::string& test_name) {
    if (condition) {
        std::cout << "[PASS] " << test_name << std::endl;
        tests_passed++;
    } else {
        std::cout << "[FAIL] " << test_name << std::endl;
        tests_failed++;
    }
}

void check_approx(double actual, double expected, double tol, const std::string& test_name) {
    double rel_error = std::abs(actual - expected) / std::max(std::abs(expected), 1e-300);
    bool pass = rel_error < tol;
    if (pass) {
        std::cout << "[PASS] " << test_name 
                  << " (actual=" << std::setprecision(15) << actual 
                  << ", expected=" << expected 
                  << ", rel_err=" << std::scientific << rel_error << ")" << std::endl;
        tests_passed++;
    } else {
        std::cout << "[FAIL] " << test_name 
                  << " (actual=" << std::setprecision(15) << actual 
                  << ", expected=" << expected 
                  << ", rel_err=" << std::scientific << rel_error << ")" << std::endl;
        tests_failed++;
    }
}

template <int W>
void test_exp_basic() {
    std::cout << "\n=== Testing exp() with MPFloat<" << W << "> ===" << std::endl;
    
    // Test exp(0) = 1
    {
        MPFloat<W> x(0.0);
        MPFloat<W> result = exp(x);
        double d = result.to_double();
        check_approx(d, 1.0, 1e-14, "exp(0) = 1");
    }
    
    // Test exp(1) ≈ e
    {
        MPFloat<W> x(1.0);
        MPFloat<W> result = exp(x);
        double d = result.to_double();
        check_approx(d, E_REF, 1e-14, "exp(1) = e");
    }
    
    // Test exp(-1) ≈ 1/e
    {
        MPFloat<W> x(-1.0);
        MPFloat<W> result = exp(x);
        double d = result.to_double();
        check_approx(d, 1.0/E_REF, 1e-14, "exp(-1) = 1/e");
    }
    
    // Test exp(2) ≈ e²
    {
        MPFloat<W> x(2.0);
        MPFloat<W> result = exp(x);
        double d = result.to_double();
        check_approx(d, E_REF * E_REF, 1e-13, "exp(2) = e^2");
    }
    
    // Test exp(0.5) ≈ sqrt(e)
    {
        MPFloat<W> x(0.5);
        MPFloat<W> result = exp(x);
        double d = result.to_double();
        check_approx(d, std::sqrt(E_REF), 1e-14, "exp(0.5) = sqrt(e)");
    }
    
    // Test exp(ln(2)) ≈ 2
    {
        MPFloat<W> x(LN2_REF);
        MPFloat<W> result = exp(x);
        double d = result.to_double();
        check_approx(d, 2.0, 1e-14, "exp(ln(2)) = 2");
    }
    
    // Test larger argument
    {
        MPFloat<W> x(10.0);
        MPFloat<W> result = exp(x);
        double d = result.to_double();
        check_approx(d, std::exp(10.0), 1e-12, "exp(10)");
    }
    
    // Test small negative argument
    {
        MPFloat<W> x(-0.001);
        MPFloat<W> result = exp(x);
        double d = result.to_double();
        check_approx(d, std::exp(-0.001), 1e-14, "exp(-0.001)");
    }
}

template <int W>
void test_log_basic() {
    std::cout << "\n=== Testing log() with MPFloat<" << W << "> ===" << std::endl;
    
    // Test log(1) = 0
    {
        MPFloat<W> x(1.0);
        MPFloat<W> result = log(x);
        double d = result.to_double();
        check_approx(d, 0.0, 1e-14, "log(1) = 0");
    }
    
    // Test log(e) = 1
    {
        MPFloat<W> x(E_REF);
        MPFloat<W> result = log(x);
        double d = result.to_double();
        check_approx(d, 1.0, 1e-13, "log(e) = 1");
    }
    
    // Test log(2) ≈ ln(2)
    {
        MPFloat<W> x(2.0);
        MPFloat<W> result = log(x);
        double d = result.to_double();
        check_approx(d, LN2_REF, 1e-14, "log(2) = ln(2)");
    }
    
    // Test log(10) ≈ ln(10)
    {
        MPFloat<W> x(10.0);
        MPFloat<W> result = log(x);
        double d = result.to_double();
        check_approx(d, LN10_REF, 1e-13, "log(10) = ln(10)");
    }
    
    // Test log(e²) = 2
    {
        MPFloat<W> x(E_REF * E_REF);
        MPFloat<W> result = log(x);
        double d = result.to_double();
        check_approx(d, 2.0, 1e-13, "log(e^2) = 2");
    }
    
    // Test log(sqrt(e)) = 0.5
    {
        MPFloat<W> x(std::sqrt(E_REF));
        MPFloat<W> result = log(x);
        double d = result.to_double();
        check_approx(d, 0.5, 1e-13, "log(sqrt(e)) = 0.5");
    }
    
    // Test log(0.5) = -ln(2)
    {
        MPFloat<W> x(0.5);
        MPFloat<W> result = log(x);
        double d = result.to_double();
        check_approx(d, -LN2_REF, 1e-14, "log(0.5) = -ln(2)");
    }
    
    // Test log(100) = 2*ln(10)
    {
        MPFloat<W> x(100.0);
        MPFloat<W> result = log(x);
        double d = result.to_double();
        check_approx(d, 2.0 * LN10_REF, 1e-13, "log(100) = 2*ln(10)");
    }
}

template <int W>
void test_roundtrip() {
    std::cout << "\n=== Testing roundtrip identities with MPFloat<" << W << "> ===" << std::endl;
    
    // Test exp(log(x)) = x for various x
    double test_values[] = {0.1, 0.5, 1.5, 2.0, 3.14159, 10.0, 100.0, 0.01};
    
    for (double val : test_values) {
        MPFloat<W> x(val);
        MPFloat<W> result = exp(log(x));
        double d = result.to_double();
        std::string name = "exp(log(" + std::to_string(val) + ")) = " + std::to_string(val);
        check_approx(d, val, 1e-12, name);
    }
    
    // Test log(exp(x)) = x for various x
    double test_values2[] = {-2.0, -1.0, -0.5, 0.0, 0.5, 1.0, 2.0, 5.0};
    
    for (double val : test_values2) {
        MPFloat<W> x(val);
        MPFloat<W> result = log(exp(x));
        double d = result.to_double();
        std::string name = "log(exp(" + std::to_string(val) + ")) = " + std::to_string(val);
        check_approx(d, val, 1e-12, name);
    }
}

template <int W>
void test_expm1_log1p() {
    std::cout << "\n=== Testing expm1() and log1p() with MPFloat<" << W << "> ===" << std::endl;
    
    // Test expm1 for small arguments
    {
        MPFloat<W> x(0.001);
        MPFloat<W> result = expm1(x);
        double d = result.to_double();
        check_approx(d, std::expm1(0.001), 1e-14, "expm1(0.001)");
    }
    
    {
        MPFloat<W> x(0.0001);
        MPFloat<W> result = expm1(x);
        double d = result.to_double();
        check_approx(d, std::expm1(0.0001), 1e-14, "expm1(0.0001)");
    }
    
    {
        MPFloat<W> x(-0.001);
        MPFloat<W> result = expm1(x);
        double d = result.to_double();
        check_approx(d, std::expm1(-0.001), 1e-14, "expm1(-0.001)");
    }
    
    // Test log1p for small arguments
    {
        MPFloat<W> x(0.001);
        MPFloat<W> result = log1p(x);
        double d = result.to_double();
        check_approx(d, std::log1p(0.001), 1e-14, "log1p(0.001)");
    }
    
    {
        MPFloat<W> x(0.0001);
        MPFloat<W> result = log1p(x);
        double d = result.to_double();
        check_approx(d, std::log1p(0.0001), 1e-14, "log1p(0.0001)");
    }
    
    {
        MPFloat<W> x(-0.001);
        MPFloat<W> result = log1p(x);
        double d = result.to_double();
        check_approx(d, std::log1p(-0.001), 1e-14, "log1p(-0.001)");
    }
}

template <int W>
void test_pow() {
    std::cout << "\n=== Testing pow() with MPFloat<" << W << "> ===" << std::endl;
    
    // Integer power tests
    {
        MPFloat<W> base(2.0);
        MPFloat<W> result = pow(base, 10);
        double d = result.to_double();
        check_approx(d, 1024.0, 1e-14, "pow(2, 10) = 1024");
    }
    
    {
        MPFloat<W> base(2.0);
        MPFloat<W> result = pow(base, -3);
        double d = result.to_double();
        check_approx(d, 0.125, 1e-14, "pow(2, -3) = 0.125");
    }
    
    {
        MPFloat<W> base(3.0);
        MPFloat<W> result = pow(base, 0);
        double d = result.to_double();
        check_approx(d, 1.0, 1e-14, "pow(3, 0) = 1");
    }
    
    // General power tests
    {
        MPFloat<W> base(2.0);
        MPFloat<W> exponent(0.5);
        MPFloat<W> result = pow(base, exponent);
        double d = result.to_double();
        check_approx(d, std::sqrt(2.0), 1e-13, "pow(2, 0.5) = sqrt(2)");
    }
    
    {
        MPFloat<W> base(E_REF);
        MPFloat<W> exponent(2.0);
        MPFloat<W> result = pow(base, exponent);
        double d = result.to_double();
        check_approx(d, E_REF * E_REF, 1e-13, "pow(e, 2) = e^2");
    }
    
    {
        MPFloat<W> base(10.0);
        MPFloat<W> exponent(0.5);
        MPFloat<W> result = pow(base, exponent);
        double d = result.to_double();
        check_approx(d, std::sqrt(10.0), 1e-13, "pow(10, 0.5) = sqrt(10)");
    }
}

template <int W>
void test_log2_log10() {
    std::cout << "\n=== Testing log2() and log10() with MPFloat<" << W << "> ===" << std::endl;
    
    // log2 tests
    {
        MPFloat<W> x(2.0);
        MPFloat<W> result = log2(x);
        double d = result.to_double();
        check_approx(d, 1.0, 1e-13, "log2(2) = 1");
    }
    
    {
        MPFloat<W> x(8.0);
        MPFloat<W> result = log2(x);
        double d = result.to_double();
        check_approx(d, 3.0, 1e-13, "log2(8) = 3");
    }
    
    {
        MPFloat<W> x(1024.0);
        MPFloat<W> result = log2(x);
        double d = result.to_double();
        check_approx(d, 10.0, 1e-13, "log2(1024) = 10");
    }
    
    // log10 tests
    {
        MPFloat<W> x(10.0);
        MPFloat<W> result = log10(x);
        double d = result.to_double();
        check_approx(d, 1.0, 1e-13, "log10(10) = 1");
    }
    
    {
        MPFloat<W> x(100.0);
        MPFloat<W> result = log10(x);
        double d = result.to_double();
        check_approx(d, 2.0, 1e-13, "log10(100) = 2");
    }
    
    {
        MPFloat<W> x(1000.0);
        MPFloat<W> result = log10(x);
        double d = result.to_double();
        check_approx(d, 3.0, 1e-13, "log10(1000) = 3");
    }
}

template <int W>
void test_edge_cases() {
    std::cout << "\n=== Testing edge cases with MPFloat<" << W << "> ===" << std::endl;
    
    // exp of very small number
    {
        MPFloat<W> x(1e-10);
        MPFloat<W> result = exp(x);
        double d = result.to_double();
        check_approx(d, std::exp(1e-10), 1e-14, "exp(1e-10)");
    }
    
    // log of number very close to 1
    {
        MPFloat<W> x(1.0000001);
        MPFloat<W> result = log(x);
        double d = result.to_double();
        check_approx(d, std::log(1.0000001), 1e-13, "log(1.0000001)");
    }
    
    // log of very small positive number
    {
        MPFloat<W> x(1e-10);
        MPFloat<W> result = log(x);
        double d = result.to_double();
        check_approx(d, std::log(1e-10), 1e-12, "log(1e-10)");
    }
    
    // log of large number
    {
        MPFloat<W> x(1e10);
        MPFloat<W> result = log(x);
        double d = result.to_double();
        check_approx(d, std::log(1e10), 1e-12, "log(1e10)");
    }
}

int main() {
    std::cout << "==========================================" << std::endl;
    std::cout << "   MPFloat Exp/Log Function Tests" << std::endl;
    std::cout << "==========================================" << std::endl;
    
    // Test with different precision levels
    // MPFloat<6> gives ~100 decimal digits
    test_exp_basic<6>();
    test_log_basic<6>();
    test_roundtrip<6>();
    test_expm1_log1p<6>();
    test_pow<6>();
    test_log2_log10<6>();
    test_edge_cases<6>();
    
    // Summary
    std::cout << "\n==========================================" << std::endl;
    std::cout << "   Test Summary" << std::endl;
    std::cout << "==========================================" << std::endl;
    std::cout << "Passed: " << tests_passed << std::endl;
    std::cout << "Failed: " << tests_failed << std::endl;
    std::cout << "Total:  " << (tests_passed + tests_failed) << std::endl;
    
    return tests_failed > 0 ? 1 : 0;
}
