/**
 * @file bench_arithmetic.cpp
 * @brief Benchmarks for MPFloat arithmetic operations.
 * 
 * Uses Kokkos timer for cross-platform performance measurement.
 * Runs both serial and parallel benchmarks.
 */

#include <mpfun/mpfun.hpp>
#include <Kokkos_Core.hpp>
#include <iostream>
#include <iomanip>
#include <chrono>

using namespace mpfun;

// Benchmark configuration
constexpr int WARMUP_ITERS = 10;
constexpr int BENCH_ITERS = 100;
constexpr int PARALLEL_N = 10000;

// Timer utility
class Timer {
public:
    void start() { start_ = std::chrono::high_resolution_clock::now(); }
    void stop() { end_ = std::chrono::high_resolution_clock::now(); }
    
    double elapsed_ms() const {
        return std::chrono::duration<double, std::milli>(end_ - start_).count();
    }
    
    double elapsed_us() const {
        return std::chrono::duration<double, std::micro>(end_ - start_).count();
    }
    
private:
    std::chrono::high_resolution_clock::time_point start_, end_;
};

// Print benchmark result
void print_result(const char* name, double total_ms, int iters) {
    double per_op_us = (total_ms * 1000.0) / iters;
    double ops_per_sec = (iters * 1000.0) / total_ms;
    
    std::cout << std::setw(30) << std::left << name
              << std::setw(12) << std::right << std::fixed << std::setprecision(3) << total_ms << " ms"
              << std::setw(12) << std::right << std::fixed << std::setprecision(3) << per_op_us << " us/op"
              << std::setw(15) << std::right << std::fixed << std::setprecision(0) << ops_per_sec << " ops/s"
              << std::endl;
}

// ============================================================================
// Serial Benchmarks (will be enabled once Agent 1 delivers)
// ============================================================================

/*
template <int W>
void bench_addition_serial() {
    MPFloat<W> a(3.14159265358979323846);
    MPFloat<W> b(2.71828182845904523536);
    MPFloat<W> c;
    
    Timer timer;
    
    // Warmup
    for (int i = 0; i < WARMUP_ITERS; ++i) {
        c = a + b;
        Kokkos::fence();
    }
    
    // Benchmark
    timer.start();
    for (int i = 0; i < BENCH_ITERS; ++i) {
        c = a + b;
    }
    Kokkos::fence();
    timer.stop();
    
    char name[64];
    snprintf(name, sizeof(name), "Addition (MPFloat<%d>)", W);
    print_result(name, timer.elapsed_ms(), BENCH_ITERS);
}

template <int W>
void bench_multiplication_serial() {
    MPFloat<W> a(3.14159265358979323846);
    MPFloat<W> b(2.71828182845904523536);
    MPFloat<W> c;
    
    Timer timer;
    
    // Warmup
    for (int i = 0; i < WARMUP_ITERS; ++i) {
        c = a * b;
        Kokkos::fence();
    }
    
    // Benchmark
    timer.start();
    for (int i = 0; i < BENCH_ITERS; ++i) {
        c = a * b;
    }
    Kokkos::fence();
    timer.stop();
    
    char name[64];
    snprintf(name, sizeof(name), "Multiplication (MPFloat<%d>)", W);
    print_result(name, timer.elapsed_ms(), BENCH_ITERS);
}

template <int W>
void bench_division_serial() {
    MPFloat<W> a(3.14159265358979323846);
    MPFloat<W> b(2.71828182845904523536);
    MPFloat<W> c;
    
    Timer timer;
    
    // Warmup
    for (int i = 0; i < WARMUP_ITERS; ++i) {
        c = a / b;
        Kokkos::fence();
    }
    
    // Benchmark
    timer.start();
    for (int i = 0; i < BENCH_ITERS; ++i) {
        c = a / b;
    }
    Kokkos::fence();
    timer.stop();
    
    char name[64];
    snprintf(name, sizeof(name), "Division (MPFloat<%d>)", W);
    print_result(name, timer.elapsed_ms(), BENCH_ITERS);
}
*/

// ============================================================================
// Parallel Benchmarks (will be enabled once Agent 1 delivers)
// ============================================================================

/*
template <int W>
void bench_addition_parallel() {
    using mp_type = MPFloat<W>;
    
    Kokkos::View<mp_type*> a("a", PARALLEL_N);
    Kokkos::View<mp_type*> b("b", PARALLEL_N);
    Kokkos::View<mp_type*> c("c", PARALLEL_N);
    
    // Initialize
    Kokkos::parallel_for("init", PARALLEL_N, KOKKOS_LAMBDA(int i) {
        a(i) = mp_type(static_cast<double>(i) + 1.0);
        b(i) = mp_type(static_cast<double>(i) + 2.0);
    });
    Kokkos::fence();
    
    Timer timer;
    
    // Warmup
    for (int iter = 0; iter < WARMUP_ITERS; ++iter) {
        Kokkos::parallel_for("add_warmup", PARALLEL_N, KOKKOS_LAMBDA(int i) {
            c(i) = a(i) + b(i);
        });
        Kokkos::fence();
    }
    
    // Benchmark
    timer.start();
    for (int iter = 0; iter < BENCH_ITERS; ++iter) {
        Kokkos::parallel_for("add_bench", PARALLEL_N, KOKKOS_LAMBDA(int i) {
            c(i) = a(i) + b(i);
        });
        Kokkos::fence();
    }
    timer.stop();
    
    char name[64];
    snprintf(name, sizeof(name), "Parallel Add (MPFloat<%d>, N=%d)", W, PARALLEL_N);
    print_result(name, timer.elapsed_ms(), BENCH_ITERS * PARALLEL_N);
}
*/

// ============================================================================
// Main
// ============================================================================

int main(int argc, char* argv[]) {
    Kokkos::initialize(argc, argv);
    
    {
        std::cout << "============================================================\n";
        std::cout << "MPFloat Arithmetic Benchmarks\n";
        std::cout << "============================================================\n";
        std::cout << "\nKokkos configuration:\n";
        Kokkos::print_configuration(std::cout);
        std::cout << "\n";
        
        std::cout << "Benchmarks will be enabled once arithmetic operations are implemented.\n";
        std::cout << "Currently just verifying infrastructure compiles and runs.\n\n";
        
        // Verify basic type creation works
        MPFloat100 x;
        std::cout << "Created MPFloat100: is_zero=" << x.is_zero() 
                  << ", capacity=" << x.capacity() 
                  << ", approx_digits=" << x.approx_digits() << "\n";
        
        /*
        // Serial benchmarks at different precisions
        std::cout << "\n--- Serial Benchmarks ---\n";
        std::cout << std::setw(30) << std::left << "Operation"
                  << std::setw(12) << std::right << "Total"
                  << std::setw(12) << std::right << "Per Op"
                  << std::setw(15) << std::right << "Throughput"
                  << std::endl;
        std::cout << std::string(69, '-') << std::endl;
        
        bench_addition_serial<6>();    // ~100 digits
        bench_addition_serial<28>();   // ~500 digits
        bench_addition_serial<56>();   // ~1000 digits
        
        bench_multiplication_serial<6>();
        bench_multiplication_serial<28>();
        bench_multiplication_serial<56>();
        
        bench_division_serial<6>();
        bench_division_serial<28>();
        bench_division_serial<56>();
        
        // Parallel benchmarks
        std::cout << "\n--- Parallel Benchmarks ---\n";
        std::cout << std::setw(30) << std::left << "Operation"
                  << std::setw(12) << std::right << "Total"
                  << std::setw(12) << std::right << "Per Op"
                  << std::setw(15) << std::right << "Throughput"
                  << std::endl;
        std::cout << std::string(69, '-') << std::endl;
        
        bench_addition_parallel<6>();
        bench_addition_parallel<28>();
        */
        
        std::cout << "\n============================================================\n";
    }
    
    Kokkos::finalize();
    return 0;
}
