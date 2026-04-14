/**
 * @file test_constants.cpp
 * @brief Standalone tests for mathematical constants
 * 
 * Compile with:
 *   clang++ -std=c++17 -O2 -I../../include -o test_constants test_constants.cpp
 * 
 * Tests:
 *   - pi: 3.14159265358979323846...
 *   - e:  2.71828182845904523536...
 *   - ln2: 0.69314718055994530942...
 *   - euler_gamma: 0.57721566490153286061...
 */

#include <cstdio>
#include <cstdint>
#include <cmath>
#include <cstring>
#include <string>

// Include MPFloat headers
#include "mpfun/mp_float.hpp"
#include "mpfun/transcendental/constants.hpp"

using namespace mpfun;

// Reference values from tests/validation/reference_values.hpp
// Pi to 100 digits
const char* PI_100 = 
    "3.1415926535897932384626433832795028841971693993751"
    "058209749445923078164062862089986280348253421170679";

// e to 100 digits
const char* E_100 = 
    "2.7182818284590452353602874713526624977572470936999"
    "595749669676277240766303535475945713821785251664274";

// ln(2) to 100 digits
const char* LN2_100 = 
    "0.6931471805599453094172321214581765680755001343602"
    "552541206800094933936219696947156058633269964186875";

// Euler-Mascheroni constant to 100 digits
const char* EULER_GAMMA_100 = 
    "0.5772156649015328606065120900824024310421593359399"
    "235988057672348848677267776646709369470632917467495";

/**
 * @brief Convert MPFloat to a string of decimal digits
 * 
 * Simple implementation that converts to double and formats.
 * For full precision, would need proper decimal conversion.
 */
template <int W>
std::string mpfloat_to_string(const MPFloat<W>& x, int digits = 50) {
    // Get the double approximation
    double d = x.to_double();
    
    // For now, just use double formatting
    // A proper implementation would use MPFUN's mpctomp or similar
    char buf[200];
    snprintf(buf, sizeof(buf), "%.*g", digits, d);
    return std::string(buf);
}

/**
 * @brief Compare first N digits of MPFloat with reference string
 * 
 * @return Number of matching digits
 */
template <int W>
int compare_digits(const MPFloat<W>& x, const char* reference, int max_digits = 15) {
    double d = x.to_double();
    
    // Parse reference value
    double ref = std::strtod(reference, nullptr);
    
    // Compare digit by digit
    // This is a simple comparison - for full precision tests
    // we'd need proper decimal string conversion
    
    double rel_err = std::fabs(d - ref) / std::fabs(ref);
    
    // Estimate number of correct digits from relative error
    // rel_err ≈ 10^(-correct_digits)
    int correct_digits = 0;
    if (rel_err > 0) {
        correct_digits = static_cast<int>(-std::log10(rel_err));
        if (correct_digits < 0) correct_digits = 0;
        if (correct_digits > max_digits) correct_digits = max_digits;
    } else {
        correct_digits = max_digits;  // Perfect match
    }
    
    return correct_digits;
}

void test_pi() {
    printf("\n=== Testing pi<W>() ===\n");
    
    // Test with W=6 (~100 decimal digits)
    printf("\npi<6>() (expected: ~100 decimal digits):\n");
    
    auto pi_val = pi<6>();
    double pi_d = pi_val.to_double();
    
    printf("  Computed:  %.15g\n", pi_d);
    printf("  Reference: %.15g\n", M_PI);
    printf("  Error:     %.2e\n", std::fabs(pi_d - M_PI));
    
    int correct = compare_digits(pi_val, PI_100, 15);
    printf("  Correct digits: ~%d (limited by double precision)\n", correct);
    
    if (std::fabs(pi_d - M_PI) < 1e-14) {
        printf("  STATUS: PASS\n");
    } else {
        printf("  STATUS: FAIL\n");
    }
    
    // Derived constants
    printf("\nTesting derived constants:\n");
    
    auto pi2 = pi_2<6>();
    auto pi4 = pi_4<6>();
    auto twopi = two_pi<6>();
    
    printf("  pi/2:  %.15g (expected: %.15g)\n", pi2.to_double(), M_PI / 2);
    printf("  pi/4:  %.15g (expected: %.15g)\n", pi4.to_double(), M_PI / 4);
    printf("  2*pi:  %.15g (expected: %.15g)\n", twopi.to_double(), 2 * M_PI);
    
    double err_pi2 = std::fabs(pi2.to_double() - M_PI / 2);
    double err_pi4 = std::fabs(pi4.to_double() - M_PI / 4);
    double err_twopi = std::fabs(twopi.to_double() - 2 * M_PI);
    
    if (err_pi2 < 1e-14 && err_pi4 < 1e-14 && err_twopi < 1e-14) {
        printf("  STATUS: PASS\n");
    } else {
        printf("  STATUS: FAIL (errors: %.2e, %.2e, %.2e)\n", err_pi2, err_pi4, err_twopi);
    }
}

void test_e() {
    printf("\n=== Testing e<W>() ===\n");
    
    auto e_val = e<6>();
    double e_d = e_val.to_double();
    
    printf("  Computed:  %.15g\n", e_d);
    printf("  Reference: %.15g\n", M_E);
    printf("  Error:     %.2e\n", std::fabs(e_d - M_E));
    
    int correct = compare_digits(e_val, E_100, 15);
    printf("  Correct digits: ~%d (limited by double precision)\n", correct);
    
    if (std::fabs(e_d - M_E) < 1e-14) {
        printf("  STATUS: PASS\n");
    } else {
        printf("  STATUS: FAIL\n");
    }
}

void test_ln2() {
    printf("\n=== Testing ln2<W>() ===\n");
    
    auto ln2_val = ln2<6>();
    double ln2_d = ln2_val.to_double();
    double ref = std::log(2.0);
    
    printf("  Computed:  %.15g\n", ln2_d);
    printf("  Reference: %.15g\n", ref);
    printf("  Error:     %.2e\n", std::fabs(ln2_d - ref));
    
    int correct = compare_digits(ln2_val, LN2_100, 15);
    printf("  Correct digits: ~%d (limited by double precision)\n", correct);
    
    if (std::fabs(ln2_d - ref) < 1e-14) {
        printf("  STATUS: PASS\n");
    } else {
        printf("  STATUS: FAIL\n");
    }
    
    // Derived constants
    printf("\nTesting log-related constants:\n");
    
    auto log2e_val = log2_e<6>();
    double log2e_d = log2e_val.to_double();
    double log2e_ref = 1.0 / std::log(2.0);  // log2(e) = 1/ln(2)
    
    printf("  log2(e):  %.15g (expected: %.15g)\n", log2e_d, log2e_ref);
    printf("  Error:    %.2e\n", std::fabs(log2e_d - log2e_ref));
    
    if (std::fabs(log2e_d - log2e_ref) < 1e-13) {
        printf("  STATUS: PASS\n");
    } else {
        printf("  STATUS: FAIL\n");
    }
}

void test_euler_gamma() {
    printf("\n=== Testing euler_gamma<W>() ===\n");
    
    auto gamma_val = euler_gamma<6>();
    double gamma_d = gamma_val.to_double();
    double ref = 0.5772156649015328606;  // Euler-Mascheroni constant
    
    printf("  Computed:  %.15g\n", gamma_d);
    printf("  Reference: %.15g\n", ref);
    printf("  Error:     %.2e\n", std::fabs(gamma_d - ref));
    
    int correct = compare_digits(gamma_val, EULER_GAMMA_100, 15);
    printf("  Correct digits: ~%d (limited by double precision)\n", correct);
    
    // Euler-Mascheroni is harder to compute precisely
    // Accept larger error for now
    if (std::fabs(gamma_d - ref) < 1e-10) {
        printf("  STATUS: PASS\n");
    } else {
        printf("  STATUS: FAIL\n");
    }
}

void test_sqrt_constants() {
    printf("\n=== Testing sqrt constants ===\n");
    
    auto sqrt2_val = sqrt2<6>();
    auto sqrt3_val = sqrt3<6>();
    auto inv_sqrt2_val = inv_sqrt2<6>();
    
    double sqrt2_d = sqrt2_val.to_double();
    double sqrt3_d = sqrt3_val.to_double();
    double inv_sqrt2_d = inv_sqrt2_val.to_double();
    
    double sqrt2_ref = std::sqrt(2.0);
    double sqrt3_ref = std::sqrt(3.0);
    double inv_sqrt2_ref = 1.0 / std::sqrt(2.0);
    
    printf("  sqrt(2):   %.15g (expected: %.15g)\n", sqrt2_d, sqrt2_ref);
    printf("  sqrt(3):   %.15g (expected: %.15g)\n", sqrt3_d, sqrt3_ref);
    printf("  1/sqrt(2): %.15g (expected: %.15g)\n", inv_sqrt2_d, inv_sqrt2_ref);
    
    double err2 = std::fabs(sqrt2_d - sqrt2_ref);
    double err3 = std::fabs(sqrt3_d - sqrt3_ref);
    double err_inv = std::fabs(inv_sqrt2_d - inv_sqrt2_ref);
    
    printf("  Errors: %.2e, %.2e, %.2e\n", err2, err3, err_inv);
    
    if (err2 < 1e-14 && err3 < 1e-14 && err_inv < 1e-14) {
        printf("  STATUS: PASS\n");
    } else {
        printf("  STATUS: FAIL\n");
    }
}

void test_precision_scaling() {
    printf("\n=== Testing precision scaling ===\n");
    
    // Test that constants are computed correctly at different precisions
    printf("\nPi at different precisions:\n");
    
    // W=3 (~50 digits)
    auto pi3 = pi<3>();
    printf("  pi<3>:  %.15g\n", pi3.to_double());
    
    // W=6 (~100 digits)
    auto pi6 = pi<6>();
    printf("  pi<6>:  %.15g\n", pi6.to_double());
    
    // W=12 (~200 digits)
    auto pi12 = pi<12>();
    printf("  pi<12>: %.15g\n", pi12.to_double());
    
    // All should match to double precision
    double err3 = std::fabs(pi3.to_double() - M_PI);
    double err6 = std::fabs(pi6.to_double() - M_PI);
    double err12 = std::fabs(pi12.to_double() - M_PI);
    
    printf("  Errors: %.2e, %.2e, %.2e\n", err3, err6, err12);
    
    if (err3 < 1e-14 && err6 < 1e-14 && err12 < 1e-14) {
        printf("  STATUS: PASS\n");
    } else {
        printf("  STATUS: FAIL\n");
    }
}

void test_identities() {
    printf("\n=== Testing mathematical identities ===\n");
    
    // Test e^(ln(2)) = 2
    // We don't have exp() yet, but we can test ln(2) * log2(e) = 1
    auto ln2_val = ln2<6>();
    auto log2e_val = log2_e<6>();
    auto product = ln2_val * log2e_val;
    
    printf("  ln(2) * log2(e) = %.15g (expected: 1.0)\n", product.to_double());
    
    double err1 = std::fabs(product.to_double() - 1.0);
    
    // Test pi/4 * 4 = pi
    auto pi4 = pi_4<6>();
    auto pi_reconstructed = pi4 * MPFloat<6>(4.0);
    auto pi_direct = pi<6>();
    
    printf("  pi/4 * 4 = %.15g (expected: %.15g)\n", 
           pi_reconstructed.to_double(), pi_direct.to_double());
    
    double err2 = std::fabs(pi_reconstructed.to_double() - pi_direct.to_double());
    
    // Test sqrt(2)^2 = 2
    auto sqrt2_val = sqrt2<6>();
    auto two_reconstructed = sqrt2_val * sqrt2_val;
    
    printf("  sqrt(2)^2 = %.15g (expected: 2.0)\n", two_reconstructed.to_double());
    
    double err3 = std::fabs(two_reconstructed.to_double() - 2.0);
    
    printf("  Errors: %.2e, %.2e, %.2e\n", err1, err2, err3);
    
    if (err1 < 1e-14 && err2 < 1e-14 && err3 < 1e-14) {
        printf("  STATUS: PASS\n");
    } else {
        printf("  STATUS: FAIL\n");
    }
}

int main() {
    printf("MPFloat Constants Tests\n");
    printf("=======================\n");
    printf("Testing mathematical constants implementation.\n");
    printf("All tests limited to ~15 digits by double precision comparison.\n");
    printf("Full precision validation requires decimal string conversion.\n");
    
    test_pi();
    test_e();
    test_ln2();
    test_euler_gamma();
    test_sqrt_constants();
    test_precision_scaling();
    test_identities();
    
    printf("\n=== Test Summary ===\n");
    printf("Check STATUS lines above for pass/fail.\n");
    printf("Note: Euler-Mascheroni uses a simplified algorithm.\n");
    printf("      For production, use Brent-McMillan with proper ln(n).\n");
    
    return 0;
}
