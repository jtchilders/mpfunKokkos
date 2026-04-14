/**
 * @file test_div.cpp
 * @brief Standalone debug test for MPFloat division
 * 
 * Compile with:
 *   clang++ -std=c++17 -I../../include -o test_div test_div.cpp
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

using namespace mpfun::detail;

constexpr int WORDS = 32;

// Helper to print an MPFloat value
void print_mpfloat(const char* name, const int64_t* a) {
    int na = extract_length(a[IDX_SIGN_LENGTH]);
    int sa = extract_sign(a[IDX_SIGN_LENGTH]);
    int64_t exp = a[IDX_EXPONENT];
    
    printf("%s: sign=%d, length=%d, exponent=%lld\n", name, sa, na, (long long)exp);
    printf("  mantissa: [");
    for (int i = 0; i < (na < 8 ? na : 8); ++i) {
        if (i > 0) printf(", ");
        printf("0x%015llx", (unsigned long long)a[IDX_MANTISSA_START + i]);
    }
    if (na > 8) printf(", ...");
    printf("]\n");
    
    // Also compute approximate double value
    double approx = 0.0;
    double mult = pow(2.0, BITS_PER_WORD * exp);
    for (int i = 0; i < (na < 4 ? na : 4); ++i) {
        approx += static_cast<double>(a[IDX_MANTISSA_START + i]) * mult;
        mult /= pow(2.0, BITS_PER_WORD);
    }
    approx *= sa;
    printf("  approx double: %.15g\n", approx);
}

// Initialize an MPFloat array
void init_mpfloat(int64_t* a, int mpnw) {
    memset(a, 0, sizeof(int64_t) * (mpnw + 10));
    a[IDX_ALLOCATED] = mpnw + 6;
    a[IDX_PRECISION] = mpnw;
}

// Convert double to MPFloat (using mpdmc pattern)
void from_double(double d, int64_t* a, int mpnw) {
    mpdmc<WORDS>(d, 0, a, mpnw);
}

// Convert MPFloat back to double (for verification)
double to_double(const int64_t* a, int mpnw) {
    double b;
    int n;
    mpmdc<WORDS>(a, b, n, mpnw);
    return b * pow(2.0, n);
}

void test_division(double num, double denom, const char* label) {
    constexpr int mpnw = 8;  // Working precision
    
    int64_t a[WORDS + 10];
    int64_t b[WORDS + 10];
    int64_t c[WORDS + 10];
    
    init_mpfloat(a, mpnw);
    init_mpfloat(b, mpnw);
    init_mpfloat(c, mpnw);
    
    printf("\n=== Test: %s ===\n", label);
    printf("Computing: %.15g / %.15g = %.15g (expected)\n", num, denom, num / denom);
    
    from_double(num, a, mpnw);
    from_double(denom, b, mpnw);
    
    printf("\nInputs:\n");
    print_mpfloat("a (numerator)", a);
    print_mpfloat("b (denominator)", b);
    
    // Perform division
    div<WORDS>(a, b, c, mpnw);
    
    printf("\nResult:\n");
    print_mpfloat("c (a/b)", c);
    
    double result = to_double(c, mpnw);
    double expected = num / denom;
    double rel_error = fabs(result - expected) / fabs(expected);
    
    printf("\nVerification:\n");
    printf("  Expected:  %.15g\n", expected);
    printf("  Got:       %.15g\n", result);
    printf("  Rel Error: %.2e\n", rel_error);
    
    if (rel_error < 1e-14) {
        printf("  STATUS: PASS\n");
    } else {
        printf("  STATUS: FAIL ***\n");
    }
}

void test_multiplication(double a_val, double b_val, const char* label) {
    constexpr int mpnw = 8;
    
    int64_t a[WORDS + 10];
    int64_t b[WORDS + 10];
    int64_t c[WORDS + 10];
    
    init_mpfloat(a, mpnw);
    init_mpfloat(b, mpnw);
    init_mpfloat(c, mpnw);
    
    printf("\n=== Mul Test: %s ===\n", label);
    printf("Computing: %.15g * %.15g = %.15g (expected)\n", a_val, b_val, a_val * b_val);
    
    from_double(a_val, a, mpnw);
    from_double(b_val, b, mpnw);
    
    mul<WORDS>(a, b, c, mpnw);
    
    double result = to_double(c, mpnw);
    double expected = a_val * b_val;
    double rel_error = (expected != 0) ? fabs(result - expected) / fabs(expected) : fabs(result);
    
    printf("  Expected:  %.15g\n", expected);
    printf("  Got:       %.15g\n", result);
    printf("  Rel Error: %.2e\n", rel_error);
    printf("  STATUS: %s\n", (rel_error < 1e-14) ? "PASS" : "FAIL ***");
}

void debug_multiply(double a_val, double b_val) {
    constexpr int mpnw = 8;
    
    int64_t a[WORDS + 10];
    int64_t b[WORDS + 10];
    int64_t c[WORDS + 10];
    
    init_mpfloat(a, mpnw);
    init_mpfloat(b, mpnw);
    init_mpfloat(c, mpnw);
    
    from_double(a_val, a, mpnw);
    from_double(b_val, b, mpnw);
    
    printf("a = %.15g:\n", a_val);
    print_mpfloat("  ", a);
    printf("b = %.15g:\n", b_val);
    print_mpfloat("  ", b);
    
    mul<WORDS>(a, b, c, mpnw);
    
    printf("a * b (expected %.15g):\n", a_val * b_val);
    print_mpfloat("  ", c);
    printf("  got double: %.15g\n", to_double(c, mpnw));
    
    // Now subtract 1
    int64_t one[WORDS + 10];
    int64_t diff[WORDS + 10];
    init_mpfloat(one, mpnw);
    init_mpfloat(diff, mpnw);
    from_double(1.0, one, mpnw);
    
    printf("one = 1.0:\n");
    print_mpfloat("  ", one);
    
    sub<WORDS>(one, c, diff, mpnw);
    printf("1 - (a*b) (expected ~0):\n");
    print_mpfloat("  ", diff);
    printf("  got double: %.15g\n", to_double(diff, mpnw));
}

void debug_newton_raphson(double denom) {
    constexpr int mpnw = 8;
    
    int64_t b[WORDS + 10];
    int64_t s0[WORDS + 10];
    int64_t s1[WORDS + 10];
    int64_t s2[WORDS + 10];
    int64_t s3[WORDS + 10];
    
    for (int i = 0; i < WORDS + 10; ++i) {
        b[i] = s0[i] = s1[i] = s2[i] = s3[i] = 0;
    }
    
    init_mpfloat(b, mpnw);
    init_mpfloat(s0, mpnw);
    init_mpfloat(s1, mpnw);
    init_mpfloat(s2, mpnw);
    init_mpfloat(s3, mpnw);
    
    s0[IDX_ALLOCATED] = mpnw + 7;
    s1[IDX_ALLOCATED] = mpnw + 7;
    s2[IDX_ALLOCATED] = mpnw + 7;
    s3[IDX_ALLOCATED] = mpnw + 7;
    
    printf("\n=== Debug Newton-Raphson for 1/%.15g ===\n", denom);
    
    from_double(denom, b, mpnw);
    printf("b = ");
    print_mpfloat("b", b);
    
    // Get initial approximation
    double db;
    int n_exp;
    mpmdc<WORDS>(b, db, n_exp, mpnw);
    printf("\nmpmdc result: db=%.15g, n_exp=%d\n", db, n_exp);
    
    double t2 = 1.0 / db;
    printf("Initial approx 1/db = %.15g, using exponent = %d\n", t2, -n_exp);
    
    mpdmc<WORDS>(t2, -n_exp, s2, mpnw);
    printf("\nInitial X_0 (approx 1/b):\n");
    print_mpfloat("s2", s2);
    printf("  approx value: %.15g\n", to_double(s2, mpnw));
    printf("  expected 1/b: %.15g\n", 1.0 / denom);
    
    // Initialize s3 = 1.0
    mpdmc<WORDS>(1.0, 0, s3, mpnw);
    printf("\ns3 = 1.0:\n");
    print_mpfloat("s3", s3);
    
    // Compute number of iterations
    double t1 = static_cast<double>(mpnw);
    int mq = static_cast<int>(LOG2_E * log(t1) + 1.0 - COMPARE_FUZZ);
    printf("\nNumber of iterations mq = %d\n", mq);
    
    int mpnw1 = 5;
    int nw1, nw2;
    
    for (int k = 1; k <= mq - 1; ++k) {
        if (k > 2) {
            mpnw1 = (2 * mpnw1 - 2 < mpnw) ? 2 * mpnw1 - 2 : mpnw;
            mpnw1 += 1;
        }
        nw1 = mpnw1;
        nw2 = mpnw1;
        
        printf("\n--- Iteration %d (mpnw1=%d) ---\n", k, mpnw1);
        
        mul<WORDS>(b, s2, s1, nw2);      // s1 = B * X_k
        printf("After mul(b, s2, s1): ");
        print_mpfloat("s1", s1);
        
        sub<WORDS>(s3, s1, s0, nw2);     // s0 = 1 - B * X_k
        printf("After sub(s3, s1, s0): ");
        print_mpfloat("s0", s0);
        
        mul<WORDS>(s2, s0, s1, nw1);     // s1 = X_k * (1 - B * X_k)
        printf("After mul(s2, s0, s1): ");
        print_mpfloat("s1", s1);
        
        add<WORDS>(s2, s1, s0, nw2);     // s0 = X_k + X_k * (1 - B * X_k)
        printf("After add(s2, s1, s0): ");
        print_mpfloat("s0", s0);
        
        mpeq<WORDS>(s0, s2, nw2);        // s2 = s0
        printf("After mpeq(s0, s2): ");
        print_mpfloat("s2", s2);
        
        printf("Current approximation to 1/b: %.15g (expected %.15g)\n", 
               to_double(s2, mpnw), 1.0/denom);
    }
}

int main() {
    printf("MPFloat Division Debug Tests\n");
    printf("============================\n");
    printf("BITS_PER_WORD = %d\n", BITS_PER_WORD);
    printf("RADIX = 2^%d = %.2e\n", BITS_PER_WORD, (double)RADIX);
    
    // Test multiplication first (division depends on it)
    printf("\n\n=== MULTIPLICATION TESTS ===\n");
    test_multiplication(2.0, 3.0, "2 * 3");
    test_multiplication(1.5, 2.0, "1.5 * 2");
    test_multiplication(1e10, 1e10, "1e10 * 1e10");
    test_multiplication(1e-10, 1e10, "1e-10 * 1e10");
    
    // Test basic divisions
    printf("\n\n=== DIVISION TESTS ===\n");
    test_division(1.0, 2.0, "1/2");
    test_division(4.0, 2.0, "4/2");
    test_division(10.0, 3.0, "10/3");
    test_division(1.0, 3.0, "1/3");
    test_division(22.0, 7.0, "22/7 (approx pi)");
    test_division(1e10, 1e5, "1e10/1e5");
    test_division(1e-10, 1e-5, "1e-10/1e-5");
    test_division(1e100, 1e50, "1e100/1e50");
    test_division(1e-100, 1e-50, "1e-100/1e-50");
    test_division(1.0, 7.0, "1/7");
    test_division(1.0, 9.0, "1/9");
    test_division(355.0, 113.0, "355/113 (better pi approx)");
    test_division(123456789.0, 987654321.0, "123456789/987654321");
    
    // Debug Newton-Raphson for a simple case
    debug_newton_raphson(2.0);
    
    // Debug the failing case
    printf("\n\n=== DEBUG FAILING CASE 1e10/1e5 ===\n");
    debug_newton_raphson(1e5);
    
    // Debug specific multiplication
    printf("\n\n=== DEBUG: 100000 * 1e-5 ===\n");
    debug_multiply(100000.0, 1e-5);
    
    return 0;
}
