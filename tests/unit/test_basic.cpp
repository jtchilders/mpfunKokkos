/**
 * @file test_basic.cpp
 * @brief Basic compilation and arithmetic tests for MPFloat type.
 * 
 * Tests verify that the library compiles and basic arithmetic works.
 */

#include <gtest/gtest.h>
#include <mpfun/mpfun.hpp>
#include <cmath>

namespace {

using namespace mpfun;

// ============================================================================
// Type Properties Tests
// ============================================================================

TEST(MPFloatBasic, TypeProperties) {
    static_assert(MPFloat100::capacity() == 6, "MPFloat100 should have 6 words");
    static_assert(MPFloat50::capacity() == 3, "MPFloat50 should have 3 words");
    static_assert(MPFloat1000::capacity() == 56, "MPFloat1000 should have 56 words");
}

TEST(MPFloatBasic, ApproxDigits) {
    EXPECT_GE(MPFloat100::approx_digits(), 100);
    EXPECT_LE(MPFloat100::approx_digits(), 120);
    EXPECT_GE(MPFloat50::approx_digits(), 50);
    EXPECT_GE(MPFloat1000::approx_digits(), 1000);
}

// ============================================================================
// Constructor Tests
// ============================================================================

TEST(MPFloatBasic, DefaultConstructor) {
    MPFloat100 x;
    EXPECT_TRUE(x.is_zero());
    EXPECT_FALSE(x.is_nan());
    EXPECT_EQ(x.sign(), 0);
}

TEST(MPFloatBasic, DoubleConstructor) {
    MPFloat100 a(3.14159265358979323846);
    double back = a.to_double();
    EXPECT_NEAR(back, 3.14159265358979323846, 1e-14);
}

TEST(MPFloatBasic, IntConstructor) {
    MPFloat100 a(42);
    EXPECT_NEAR(a.to_double(), 42.0, 1e-10);
}

TEST(MPFloatBasic, NegativeConstruction) {
    MPFloat100 a(-123.456);
    EXPECT_TRUE(a.is_negative());
    EXPECT_NEAR(a.to_double(), -123.456, 1e-10);
}

TEST(MPFloatBasic, ZeroState) {
    MPFloat<6> x;
    EXPECT_TRUE(x.is_zero());
    EXPECT_FALSE(x.is_negative());
    EXPECT_FALSE(x.is_positive());
}

// ============================================================================
// Accessor Tests
// ============================================================================

TEST(MPFloatBasic, DataAccess) {
    MPFloat100 x;
    auto* data = x.data();
    EXPECT_NE(data, nullptr);
    EXPECT_EQ(data[0], 6);  // capacity
    EXPECT_EQ(data[1], 6);  // working precision
    EXPECT_EQ(data[2], 0);  // zero
}

TEST(MPFloatBasic, SetNaN) {
    MPFloat100 x;
    EXPECT_FALSE(x.is_nan());
    x.set_nan();
    EXPECT_TRUE(x.is_nan());
}

TEST(MPFloatBasic, SignAccessor) {
    MPFloat100 pos(42.0);
    MPFloat100 neg(-42.0);
    MPFloat100 zero;
    
    EXPECT_EQ(pos.sign(), 1);
    EXPECT_EQ(neg.sign(), -1);
    EXPECT_EQ(zero.sign(), 0);
}

// ============================================================================
// Negation Tests
// ============================================================================

TEST(MPFloatBasic, Negation) {
    MPFloat100 a(123.456);
    MPFloat100 b = -a;
    EXPECT_NEAR(b.to_double(), -123.456, 1e-10);
    
    MPFloat100 c = -b;
    EXPECT_NEAR(c.to_double(), 123.456, 1e-10);
}

// ============================================================================
// Addition Tests
// ============================================================================

TEST(MPFloatArithmetic, AddPositives) {
    MPFloat100 a(100.0);
    MPFloat100 b(23.5);
    MPFloat100 c = a + b;
    EXPECT_NEAR(c.to_double(), 123.5, 1e-10);
}

TEST(MPFloatArithmetic, AddNegatives) {
    MPFloat100 a(-100.0);
    MPFloat100 b(-23.5);
    MPFloat100 c = a + b;
    EXPECT_NEAR(c.to_double(), -123.5, 1e-10);
}

TEST(MPFloatArithmetic, AddMixed) {
    MPFloat100 a(100.0);
    MPFloat100 b(-23.5);
    MPFloat100 c = a + b;
    EXPECT_NEAR(c.to_double(), 76.5, 1e-10);
}

TEST(MPFloatArithmetic, AddZero) {
    MPFloat100 a(42.0);
    MPFloat100 zero;
    MPFloat100 c = a + zero;
    EXPECT_NEAR(c.to_double(), 42.0, 1e-10);
}

// ============================================================================
// Subtraction Tests
// ============================================================================

TEST(MPFloatArithmetic, SubPositives) {
    MPFloat100 a(100.0);
    MPFloat100 b(23.5);
    MPFloat100 c = a - b;
    EXPECT_NEAR(c.to_double(), 76.5, 1e-10);
}

TEST(MPFloatArithmetic, SubToZero) {
    MPFloat100 a(42.0);
    MPFloat100 b(42.0);
    MPFloat100 c = a - b;
    EXPECT_TRUE(c.is_zero());
}

// ============================================================================
// Multiplication Tests
// ============================================================================

TEST(MPFloatArithmetic, MulPositives) {
    MPFloat100 a(12.0);
    MPFloat100 b(3.5);
    MPFloat100 c = a * b;
    EXPECT_NEAR(c.to_double(), 42.0, 1e-10);
}

TEST(MPFloatArithmetic, MulByZero) {
    MPFloat100 a(42.0);
    MPFloat100 zero;
    MPFloat100 c = a * zero;
    EXPECT_TRUE(c.is_zero());
}

TEST(MPFloatArithmetic, MulNegatives) {
    MPFloat100 a(-6.0);
    MPFloat100 b(-7.0);
    MPFloat100 c = a * b;
    EXPECT_NEAR(c.to_double(), 42.0, 1e-10);
}

TEST(MPFloatArithmetic, MulMixed) {
    MPFloat100 a(6.0);
    MPFloat100 b(-7.0);
    MPFloat100 c = a * b;
    EXPECT_NEAR(c.to_double(), -42.0, 1e-10);
}

// ============================================================================
// Division Tests
// ============================================================================

TEST(MPFloatArithmetic, DivExact) {
    MPFloat100 a(100.0);
    MPFloat100 b(4.0);
    MPFloat100 c = a / b;
    EXPECT_NEAR(c.to_double(), 25.0, 1e-10);
}

TEST(MPFloatArithmetic, DivNonExact) {
    MPFloat100 a(1.0);
    MPFloat100 b(3.0);
    MPFloat100 c = a / b;
    EXPECT_NEAR(c.to_double(), 1.0/3.0, 1e-14);
}

TEST(MPFloatArithmetic, DivByZero) {
    MPFloat100 a(42.0);
    MPFloat100 zero;
    MPFloat100 c = a / zero;
    EXPECT_TRUE(c.is_nan());
}

// ============================================================================
// Comparison Tests
// ============================================================================

TEST(MPFloatComparison, LessThan) {
    MPFloat100 a(10.0);
    MPFloat100 b(20.0);
    EXPECT_TRUE(a < b);
    EXPECT_FALSE(b < a);
}

TEST(MPFloatComparison, Equal) {
    MPFloat100 a(42.0);
    MPFloat100 b(42.0);
    EXPECT_TRUE(a == b);
    EXPECT_FALSE(a != b);
}

TEST(MPFloatComparison, NotEqual) {
    MPFloat100 a(10.0);
    MPFloat100 b(20.0);
    EXPECT_TRUE(a != b);
    EXPECT_FALSE(a == b);
}

TEST(MPFloatComparison, LessOrEqual) {
    MPFloat100 a(10.0);
    MPFloat100 b(20.0);
    MPFloat100 c(10.0);
    EXPECT_TRUE(a <= b);
    EXPECT_TRUE(a <= c);
    EXPECT_FALSE(b <= a);
}

TEST(MPFloatComparison, GreaterOrEqual) {
    MPFloat100 a(10.0);
    MPFloat100 b(20.0);
    MPFloat100 c(10.0);
    EXPECT_TRUE(b >= a);
    EXPECT_TRUE(a >= c);
    EXPECT_FALSE(a >= b);
}

// ============================================================================
// Large/Small Number Tests
// ============================================================================

TEST(MPFloatRange, LargeNumbers) {
    MPFloat100 a(1e50);
    MPFloat100 b(2e50);
    MPFloat100 c = a + b;
    EXPECT_NEAR(c.to_double(), 3e50, 1e40);
}

TEST(MPFloatRange, SmallNumbers) {
    MPFloat100 a(1e-50);
    MPFloat100 b(2e-50);
    MPFloat100 c = a + b;
    EXPECT_NEAR(c.to_double(), 3e-50, 1e-60);
}

// ============================================================================
// Abs Function Tests
// ============================================================================

TEST(MPFloatFunctions, Abs) {
    MPFloat100 a(-42.0);
    MPFloat100 b = abs(a);
    EXPECT_NEAR(b.to_double(), 42.0, 1e-10);
    EXPECT_TRUE(b.is_positive());
}

TEST(MPFloatFunctions, AbsOfPositive) {
    MPFloat100 a(42.0);
    MPFloat100 b = abs(a);
    EXPECT_NEAR(b.to_double(), 42.0, 1e-10);
}

// ============================================================================
// Compound Assignment Tests
// ============================================================================

TEST(MPFloatArithmetic, CompoundAdd) {
    MPFloat100 a(10.0);
    a += MPFloat100(5.0);
    EXPECT_NEAR(a.to_double(), 15.0, 1e-10);
}

TEST(MPFloatArithmetic, CompoundSub) {
    MPFloat100 a(10.0);
    a -= MPFloat100(3.0);
    EXPECT_NEAR(a.to_double(), 7.0, 1e-10);
}

TEST(MPFloatArithmetic, CompoundMul) {
    MPFloat100 a(10.0);
    a *= MPFloat100(3.0);
    EXPECT_NEAR(a.to_double(), 30.0, 1e-10);
}

TEST(MPFloatArithmetic, CompoundDiv) {
    MPFloat100 a(30.0);
    a /= MPFloat100(3.0);
    EXPECT_NEAR(a.to_double(), 10.0, 1e-10);
}

// ============================================================================
// Kokkos Compatibility Tests
// ============================================================================

TEST(MPFloatBasic, KokkosTypeTraits) {
    using mp_type = MPFloat100;
    static_assert(std::is_default_constructible_v<mp_type>,
                  "MPFloat should be default constructible");
}

} // namespace

int main(int argc, char** argv) {
    Kokkos::initialize(argc, argv);
    ::testing::InitGoogleTest(&argc, argv);
    int result = RUN_ALL_TESTS();
    Kokkos::finalize();
    return result;
}
