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
# Standalone Tests (Phase 2: Transcendentals)
# ============================================================================
# These are self-contained tests that don't require GTest.
# They return 0 on success, non-zero on failure.

# Core arithmetic test
add_executable(test_div standalone/test_div.cpp)
target_link_libraries(test_div PRIVATE mpfun::mpfun)
target_compile_features(test_div PRIVATE cxx_std_17)
add_test(NAME standalone_div COMMAND test_div)

# Constants test
add_executable(test_constants standalone/test_constants.cpp)
target_link_libraries(test_constants PRIVATE mpfun::mpfun)
target_compile_features(test_constants PRIVATE cxx_std_17)
add_test(NAME standalone_constants COMMAND test_constants)

# Exp/Log test
add_executable(test_exp_log standalone/test_exp_log.cpp)
target_link_libraries(test_exp_log PRIVATE mpfun::mpfun)
target_compile_features(test_exp_log PRIVATE cxx_std_17)
add_test(NAME standalone_exp_log COMMAND test_exp_log)

# Trig test
add_executable(test_trig standalone/test_trig.cpp)
target_link_libraries(test_trig PRIVATE mpfun::mpfun)
target_compile_features(test_trig PRIVATE cxx_std_17)
add_test(NAME standalone_trig COMMAND test_trig)

# Hyperbolic test
add_executable(test_hyperbolic standalone/test_hyperbolic.cpp)
target_link_libraries(test_hyperbolic PRIVATE mpfun::mpfun)
target_compile_features(test_hyperbolic PRIVATE cxx_std_17)
add_test(NAME standalone_hyperbolic COMMAND test_hyperbolic)

# ============================================================================
# Validation Tests (against high-precision reference)
# ============================================================================

# Add validation subdirectory with Python-generated test data
add_subdirectory(validation)
