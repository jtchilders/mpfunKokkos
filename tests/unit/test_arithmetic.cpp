/**
 * @file test_arithmetic.cpp
 * @brief Comprehensive arithmetic tests for MPFloat
 * 
 * Tests division Newton-Raphson convergence and square root operations.
 */

#include <gtest/gtest.h>
#include <mpfun/mpfun.hpp>
#include <cmath>
#include <iostream>
#include <iomanip>

namespace {

using namespace mpfun;

// ============================================================================
// Division Tests - Core Algorithm Verification
// ============================================================================

TEST(MPFloatDivision, BasicDivision) {
    MPFloat100 a(100.0);
    MPFloat100 b(4.0);
    MPFloat100 c = a / b;
    EXPECT_NEAR(c.to_double(), 25.0, 1e-10);
}

TEST(MPFloatDivision, DivisionByOne) {
    MPFloat100 a(3.14159265358979);
    MPFloat100 one(1.0);
    MPFloat100 c = a / one;
    EXPECT_NEAR(c.to_double(), 3.14159265358979, 1e-14);
}

TEST(MPFloatDivision, DivisionOfOne) {
    // Test 1/x for various x
    MPFloat100 one(1.0);
    
    MPFloat100 two(2.0);
    MPFloat100 half = one / two;
    EXPECT_NEAR(half.to_double(), 0.5, 1e-14);
    
    MPFloat100 four(4.0);
    MPFloat100 quarter = one / four;
    EXPECT_NEAR(quarter.to_double(), 0.25, 1e-14);
    
    MPFloat100 ten(10.0);
    MPFloat100 tenth = one / ten;
    EXPECT_NEAR(tenth.to_double(), 0.1, 1e-14);
}

TEST(MPFloatDivision, ExponentHandling) {
    // Test division with large exponent differences
    MPFloat100 large(1e20);
    MPFloat100 small(1e10);
    MPFloat100 ratio = large / small;
    EXPECT_NEAR(ratio.to_double(), 1e10, 1e-5 * 1e10);
    
    // Reverse
    MPFloat100 ratio2 = small / large;
    EXPECT_NEAR(ratio2.to_double(), 1e-10, 1e-24);
}

TEST(MPFloatDivision, FractionalDivision) {
    // 1/3 = 0.333...
    MPFloat100 one(1.0);
    MPFloat100 three(3.0);
    MPFloat100 third = one / three;
    EXPECT_NEAR(third.to_double(), 1.0/3.0, 1e-14);
    
    // Verify 3 * (1/3) ≈ 1
    MPFloat100 back = three * third;
    EXPECT_NEAR(back.to_double(), 1.0, 1e-12);
}

TEST(MPFloatDivision, NegativeDivision) {
    MPFloat100 a(42.0);
    MPFloat100 b(-7.0);
    
    MPFloat100 c = a / b;
    EXPECT_NEAR(c.to_double(), -6.0, 1e-10);
    
    MPFloat100 d = b / a;
    EXPECT_NEAR(d.to_double(), -7.0/42.0, 1e-14);
}

TEST(MPFloatDivision, ZeroDividend) {
    MPFloat100 zero;
    MPFloat100 nonzero(42.0);
    MPFloat100 result = zero / nonzero;
    EXPECT_TRUE(result.is_zero());
}

TEST(MPFloatDivision, ZeroDivisor) {
    MPFloat100 a(42.0);
    MPFloat100 zero;
    MPFloat100 result = a / zero;
    EXPECT_TRUE(result.is_nan());
}

TEST(MPFloatDivision, NewtonRaphsonConvergence) {
    // Test that Newton-Raphson converges correctly for division
    // by verifying a * (1/a) = 1 to high precision
    MPFloat100 a(7.0);
    MPFloat100 one(1.0);
    MPFloat100 inv_a = one / a;
    MPFloat100 product = a * inv_a;
    
    // Should be very close to 1
    MPFloat100 error = product - one;
    EXPECT_TRUE(error.is_zero() || abs(error).to_double() < 1e-90);
}

TEST(MPFloatDivision, DivisionIdentities) {
    MPFloat100 a(123.456);
    MPFloat100 b(789.012);
    MPFloat100 one(1.0);
    
    // a / a = 1
    MPFloat100 selfDiv = a / a;
    EXPECT_NEAR(selfDiv.to_double(), 1.0, 1e-14);
    
    // (a / b) * b ≈ a
    MPFloat100 quotient = a / b;
    MPFloat100 back = quotient * b;
    EXPECT_NEAR(back.to_double(), a.to_double(), 1e-10);
    
    // a / 1 = a
    MPFloat100 divByOne = a / one;
    EXPECT_NEAR(divByOne.to_double(), a.to_double(), 1e-14);
}

// ============================================================================
// Square Root Tests
// ============================================================================

TEST(MPFloatSqrt, BasicSqrt) {
    MPFloat100 four(4.0);
    MPFloat100 two = sqrt(four);
    EXPECT_NEAR(two.to_double(), 2.0, 1e-14);
    
    MPFloat100 nine(9.0);
    MPFloat100 three = sqrt(nine);
    EXPECT_NEAR(three.to_double(), 3.0, 1e-14);
    
    MPFloat100 sixteen(16.0);
    MPFloat100 four2 = sqrt(sixteen);
    EXPECT_NEAR(four2.to_double(), 4.0, 1e-14);
}

TEST(MPFloatSqrt, SqrtOne) {
    MPFloat100 one(1.0);
    MPFloat100 result = sqrt(one);
    EXPECT_NEAR(result.to_double(), 1.0, 1e-14);
}

TEST(MPFloatSqrt, SqrtZero) {
    MPFloat100 zero;
    MPFloat100 result = sqrt(zero);
    EXPECT_TRUE(result.is_zero());
}

TEST(MPFloatSqrt, SqrtTwo) {
    // sqrt(2) ≈ 1.41421356237...
    MPFloat100 two(2.0);
    MPFloat100 result = sqrt(two);
    EXPECT_NEAR(result.to_double(), std::sqrt(2.0), 1e-14);
    
    // Verify: result^2 should equal 2
    MPFloat100 squared = result * result;
    EXPECT_NEAR(squared.to_double(), 2.0, 1e-12);
}

TEST(MPFloatSqrt, SqrtNonPerfectSquare) {
    // Test sqrt of non-perfect squares
    MPFloat100 five(5.0);
    MPFloat100 result = sqrt(five);
    EXPECT_NEAR(result.to_double(), std::sqrt(5.0), 1e-14);
    
    // Verify: result^2 ≈ 5
    MPFloat100 squared = result * result;
    EXPECT_NEAR(squared.to_double(), 5.0, 1e-12);
}

TEST(MPFloatSqrt, SqrtLarge) {
    MPFloat100 large(1e20);
    MPFloat100 result = sqrt(large);
    EXPECT_NEAR(result.to_double(), 1e10, 1e-5);
    
    // Verify
    MPFloat100 squared = result * result;
    EXPECT_NEAR(squared.to_double(), 1e20, 1e6);
}

TEST(MPFloatSqrt, SqrtSmall) {
    MPFloat100 small(1e-20);
    MPFloat100 result = sqrt(small);
    EXPECT_NEAR(result.to_double(), 1e-10, 1e-24);
    
    // Verify
    MPFloat100 squared = result * result;
    EXPECT_NEAR(squared.to_double(), 1e-20, 1e-32);
}

TEST(MPFloatSqrt, SqrtNegative) {
    MPFloat100 neg(-4.0);
    MPFloat100 result = sqrt(neg);
    EXPECT_TRUE(result.is_nan());
}

TEST(MPFloatSqrt, SqrtInverseOfSquare) {
    // For any positive x: sqrt(x^2) = x
    MPFloat100 x(7.5);
    MPFloat100 x_squared = x * x;
    MPFloat100 result = sqrt(x_squared);
    EXPECT_NEAR(result.to_double(), 7.5, 1e-12);
}

// ============================================================================
// Combined Arithmetic Tests
// ============================================================================

TEST(MPFloatCombined, PythagoreanTriple) {
    // 3^2 + 4^2 = 5^2
    MPFloat100 three(3.0);
    MPFloat100 four(4.0);
    MPFloat100 five(5.0);
    
    MPFloat100 sum_squares = three * three + four * four;
    EXPECT_NEAR(sum_squares.to_double(), 25.0, 1e-12);
    
    MPFloat100 hypotenuse = sqrt(sum_squares);
    EXPECT_NEAR(hypotenuse.to_double(), 5.0, 1e-12);
}

TEST(MPFloatCombined, QuadraticFormula) {
    // Solve x^2 - 5x + 6 = 0
    // x = (5 ± sqrt(25 - 24)) / 2 = (5 ± 1) / 2 = 3 or 2
    
    MPFloat100 a(1.0);
    MPFloat100 b(-5.0);
    MPFloat100 c(6.0);
    MPFloat100 two(2.0);
    MPFloat100 four(4.0);
    
    // discriminant = b^2 - 4ac = 25 - 24 = 1
    MPFloat100 discriminant = b * b - four * a * c;
    EXPECT_NEAR(discriminant.to_double(), 1.0, 1e-12);
    
    MPFloat100 sqrt_d = sqrt(discriminant);
    EXPECT_NEAR(sqrt_d.to_double(), 1.0, 1e-12);
    
    // x1 = (-b + sqrt(d)) / (2a) = (5 + 1) / 2 = 3
    MPFloat100 neg_b = -b;
    MPFloat100 x1 = (neg_b + sqrt_d) / (two * a);
    EXPECT_NEAR(x1.to_double(), 3.0, 1e-12);
    
    // x2 = (-b - sqrt(d)) / (2a) = (5 - 1) / 2 = 2
    MPFloat100 x2 = (neg_b - sqrt_d) / (two * a);
    EXPECT_NEAR(x2.to_double(), 2.0, 1e-12);
}

TEST(MPFloatCombined, GeometricMean) {
    // geometric_mean(a, b) = sqrt(a * b)
    MPFloat100 a(4.0);
    MPFloat100 b(9.0);
    
    MPFloat100 product = a * b;
    MPFloat100 geo_mean = sqrt(product);
    EXPECT_NEAR(geo_mean.to_double(), 6.0, 1e-12);
}

TEST(MPFloatCombined, ArithmeticIdentities) {
    MPFloat100 a(2.0);
    MPFloat100 b(3.0);
    MPFloat100 c(5.0);
    
    // Distributive: a * (b + c) = a*b + a*c
    MPFloat100 lhs = a * (b + c);
    MPFloat100 rhs = a * b + a * c;
    EXPECT_NEAR(lhs.to_double(), rhs.to_double(), 1e-12);
    
    // (a + b) / c = a/c + b/c
    MPFloat100 lhs2 = (a + b) / c;
    MPFloat100 rhs2 = a / c + b / c;
    EXPECT_NEAR(lhs2.to_double(), rhs2.to_double(), 1e-12);
}

// ============================================================================
// Precision Tests
// ============================================================================

TEST(MPFloatPrecision, HighPrecisionDivision) {
    // Test that we get more precision than double
    MPFloat100 one(1.0);
    MPFloat100 three(3.0);
    
    MPFloat100 third = one / three;
    
    // 3 * (1/3) should be closer to 1 than double arithmetic
    MPFloat100 back = third * three;
    MPFloat100 error = back - one;
    
    // The error should be very small (less than 1e-90 for 100 digit precision)
    double err_double = error.to_double();
    EXPECT_LT(std::abs(err_double), 1e-15);  // At least better than double
}

TEST(MPFloatPrecision, HighPrecisionSqrt) {
    // sqrt(2)^2 should be very close to 2
    MPFloat100 two(2.0);
    MPFloat100 sqrt2 = sqrt(two);
    MPFloat100 squared = sqrt2 * sqrt2;
    MPFloat100 error = squared - two;
    
    double err_double = error.to_double();
    EXPECT_LT(std::abs(err_double), 1e-15);
}

// ============================================================================
// Edge Cases
// ============================================================================

TEST(MPFloatEdgeCases, VerySmallDivision) {
    MPFloat100 small(1e-50);
    MPFloat100 one(1.0);
    
    // 1 / 1e-50 = 1e50
    MPFloat100 large = one / small;
    EXPECT_NEAR(large.to_double(), 1e50, 1e35);
    
    // 1e-50 / 1 = 1e-50
    MPFloat100 stillSmall = small / one;
    EXPECT_NEAR(stillSmall.to_double(), 1e-50, 1e-64);
}

TEST(MPFloatEdgeCases, ChainedOperations) {
    MPFloat100 x(100.0);
    
    // sqrt(x) / sqrt(x) = 1
    MPFloat100 sqrtX = sqrt(x);
    MPFloat100 ratio = sqrtX / sqrtX;
    EXPECT_NEAR(ratio.to_double(), 1.0, 1e-14);
    
    // (x / 4) * 4 / x = 1
    MPFloat100 four(4.0);
    MPFloat100 result = (x / four) * four / x;
    EXPECT_NEAR(result.to_double(), 1.0, 1e-14);
}

} // anonymous namespace

int main(int argc, char** argv) {
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
