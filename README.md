# Kokkos-MPFloat: Arbitrary Precision Floating Point for Heterogeneous HPC

A Kokkos-based port of David Bailey's MPFUN2020 arbitrary precision library, designed to run on CPUs, NVIDIA GPUs (CUDA), AMD GPUs (HIP), and Intel GPUs (SYCL).

## Status

🚧 **Under Development** 🚧

## Goals

- **Portable**: Single codebase for CPU and GPU via Kokkos
- **High Performance**: Optimized for massively parallel workloads
- **Arbitrary Precision**: Compile-time configurable precision (50 to 10,000+ digits)
- **GPU-Compatible**: All core functions marked `KOKKOS_INLINE_FUNCTION`

## Quick Example

```cpp
#include <mpfun/mpfun.hpp>
#include <Kokkos_Core.hpp>

using mp100 = mpfun::MPFloat<6>;  // ~100 digits

int main() {
    Kokkos::initialize();
    {
        mp100 pi = mpfun::pi<6>();
        mp100 result;
        
        Kokkos::parallel_for(1, KOKKOS_LAMBDA(int) {
            mp100 x = mpfun::sin(pi / mp100(6.0));
            // x ≈ 0.5 to 100 digits
        });
    }
    Kokkos::finalize();
}
```

## Documentation

- [Project Plan](docs/PLAN.md) - Architecture, phases, timeline
- [Agent Tasks](docs/AGENT_TASKS.md) - Development task breakdown
- [Interface Contracts](docs/INTERFACES.md) - API specifications

## Source

Based on MPFUN2020 by David H. Bailey (LBNL):
- [MPFUN2020 Software](https://www.davidhbailey.com/dhbsoftware/)
- [Technical Paper](https://www.davidhbailey.com/dhbpapers/mpfun2020.pdf)

## License

TBD (Original MPFUN2020 is BSD-like with commercial use restrictions)

## Acknowledgments

- David H. Bailey for MPFUN2020
- Kokkos team at Sandia National Labs
