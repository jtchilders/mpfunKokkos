# ============================================================================
# Benchmark Configuration
# ============================================================================

# Benchmarks using Kokkos timing utilities
# (Google Benchmark can be added later if needed)

add_executable(bench_arithmetic bench_arithmetic.cpp)
target_link_libraries(bench_arithmetic PRIVATE mpfun::mpfun)

# More benchmarks to be added as implementations are delivered
# add_executable(bench_transcendental bench_transcendental.cpp)
# target_link_libraries(bench_transcendental PRIVATE mpfun::mpfun)
