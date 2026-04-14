# ============================================================================
# Test Configuration
# ============================================================================

# Try to find GTest, fetch if not found
include(FetchContent)

find_package(GTest QUIET)
if(NOT GTest_FOUND)
    message(STATUS "GTest not found, fetching from GitHub...")
    FetchContent_Declare(
        googletest
        GIT_REPOSITORY https://github.com/google/googletest.git
        GIT_TAG v1.14.0
    )
    # For Windows: Prevent overriding the parent project's compiler/linker settings
    set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
    FetchContent_MakeAvailable(googletest)
endif()

include(GoogleTest)

# ============================================================================
# Test Utility Library
# ============================================================================

add_library(mpfun_test_utils INTERFACE)
target_include_directories(mpfun_test_utils INTERFACE 
    ${CMAKE_CURRENT_SOURCE_DIR}
)

# ============================================================================
# Unit Tests
# ============================================================================

# Basic test that the library compiles
add_executable(test_basic unit/test_basic.cpp)
target_link_libraries(test_basic PRIVATE mpfun::mpfun GTest::gtest_main)
gtest_discover_tests(test_basic)

# Comprehensive arithmetic tests (division, sqrt)
add_executable(test_arithmetic unit/test_arithmetic.cpp)
target_link_libraries(test_arithmetic PRIVATE mpfun::mpfun GTest::gtest_main)
gtest_discover_tests(test_arithmetic)

# Add more tests as needed
# add_executable(test_add unit/test_add.cpp)
# target_link_libraries(test_add PRIVATE mpfun::mpfun GTest::gtest_main)
# gtest_discover_tests(test_add)

# add_executable(test_mul unit/test_mul.cpp)
# target_link_libraries(test_mul PRIVATE mpfun::mpfun GTest::gtest_main)
# gtest_discover_tests(test_mul)

# ============================================================================
# Validation Tests (against Fortran reference)
# ============================================================================

# add_executable(test_validation validation/test_against_reference.cpp)
# target_link_libraries(test_validation PRIVATE mpfun::mpfun GTest::gtest_main)
# gtest_discover_tests(test_validation)
