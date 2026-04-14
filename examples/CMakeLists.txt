# ============================================================================
# Examples
# ============================================================================

add_executable(basic_usage basic_usage.cpp)
target_link_libraries(basic_usage PRIVATE mpfun::mpfun)

# More examples to be added
# add_executable(pi_computation pi_computation.cpp)
# target_link_libraries(pi_computation PRIVATE mpfun::mpfun)
