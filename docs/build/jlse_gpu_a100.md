# Building and Testing on JLSE gpu_a100 Nodes

## Hardware

- **Nodes:** `gpu06`, `gpu07`
- **GPU:** NVIDIA A100-PCIE-40GB (40 GB HBM2, sm_80 / Ampere80)
- **Queue:** `gpu_a100`
- **Scheduler:** Cobalt

## Prerequisites

Jobs must be submitted via Cobalt — direct SSH to compute nodes is not available without an active allocation.

## Modules

```bash
module load cmake/3.27.2 gcc/12.2.0
module swap cuda cuda/12.3.0
```

## Step 1: Build Kokkos for A100

Clone Kokkos and build from a login node. The clone must happen on the login node (internet access); the compile job runs on the compute node.

### Clone Kokkos (login node)

```bash
KOKKOS_BASE=$HOME/kokkos_a100/kokkos-4.1.00/Kokkos_ARCH_AMPERE80/Release
mkdir -p $KOKKOS_BASE
cd $KOKKOS_BASE
git clone --branch 4.1.00 --depth 1 https://github.com/kokkos/kokkos.git
```

### Build script: `build_kokkos_a100.sh`

```bash
#!/bin/bash
set -e

module load cmake/3.27.2 gcc/12.2.0
module swap cuda cuda/12.3.0

BASEPATH=$HOME/kokkos_a100/kokkos-4.1.00/Kokkos_ARCH_AMPERE80/Release
KOKKOS_SRC=$BASEPATH/kokkos
KOKKOS_BUILD=$KOKKOS_SRC/build
KOKKOS_INSTALL=$KOKKOS_SRC/install

mkdir -p $KOKKOS_BUILD
cd $KOKKOS_BUILD

cmake $KOKKOS_SRC \
    -DCMAKE_INSTALL_PREFIX=$KOKKOS_INSTALL \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_CXX_STANDARD=17 \
    -DKokkos_ENABLE_CUDA=ON \
    -DKokkos_ARCH_AMPERE80=ON \
    -DKokkos_ENABLE_SERIAL=ON \
    -DKokkos_ENABLE_OPENMP=OFF

make -j$(nproc)
make install
```

### Submit the Kokkos build job

```bash
chmod +x build_kokkos_a100.sh
qsub -A <project> -t 60 -q gpu_a100 -n 1 \
  --output $HOME/kokkos_a100/kokkos-4.1.00/Kokkos_ARCH_AMPERE80/Release/kokkos_build.out \
  --error  $HOME/kokkos_a100/kokkos-4.1.00/Kokkos_ARCH_AMPERE80/Release/kokkos_build.err \
  build_kokkos_a100.sh
```

The install lands at:
```
$HOME/kokkos_a100/kokkos-4.1.00/Kokkos_ARCH_AMPERE80/Release/kokkos/install/
```

### Environment setup script

A convenience `setup.sh` is provided at the install root to load modules and set environment variables:

```bash
source $HOME/kokkos_a100/kokkos-4.1.00/Kokkos_ARCH_AMPERE80/Release/setup.sh
```

## Step 2: Build mpfunKokkos

### Build script: `build_mpfun_a100.sh`

```bash
#!/bin/bash
set -e

module load cmake/3.27.2 gcc/12.2.0
module swap cuda cuda/12.3.0

KOKKOS_INSTALL=$HOME/kokkos_a100/kokkos-4.1.00/Kokkos_ARCH_AMPERE80/Release/kokkos/install
MPFUN_SRC=$HOME/mpfun_testarea/mpfunKokkos
BUILD_DIR=$MPFUN_SRC/build_a100_cuda

mkdir -p $BUILD_DIR
cd $BUILD_DIR

cmake $MPFUN_SRC \
    -DCMAKE_PREFIX_PATH=$KOKKOS_INSTALL \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_CXX_STANDARD=17 \
    -DMPFUN_BUILD_TESTS=ON \
    -DMPFUN_BUILD_BENCHMARKS=OFF \
    -DMPFUN_BUILD_EXAMPLES=OFF \
    -DMPFUN_ENABLE_ASSERTIONS=ON

make -j$(nproc)
```

### Submit the mpfunKokkos build job

```bash
chmod +x build_mpfun_a100.sh
qsub -A <project> -t 60 -q gpu_a100 -n 1 \
  --output $HOME/mpfun_testarea/mpfunKokkos/build_a100_cuda.out \
  --error  $HOME/mpfun_testarea/mpfunKokkos/build_a100_cuda.err \
  build_mpfun_a100.sh
```

## Step 3: Run Tests

### Run script: `run_tests_a100.sh`

```bash
#!/bin/bash
set -e

module load cmake/3.27.2 gcc/12.2.0
module swap cuda cuda/12.3.0

BUILD=$HOME/mpfun_testarea/mpfunKokkos/build_a100_cuda

cd $BUILD
ctest --output-on-failure
```

### Submit the test job

```bash
chmod +x run_tests_a100.sh
qsub -A <project> -t 30 -q gpu_a100 -n 1 \
  --output $HOME/mpfun_testarea/mpfunKokkos/test_run_a100.out \
  --error  $HOME/mpfun_testarea/mpfunKokkos/test_run_a100.err \
  run_tests_a100.sh
```

## Test Results (as of 2026-04-15)

70 ctest tests, all passing:

| Test | Count | Notes |
|------|-------|-------|
| `MPFloatBasic` | 38 | GTest: type traits, constructors, arithmetic, comparisons |
| `MPFloatArithmetic` | 27 | GTest: division, sqrt, combined, precision |
| `standalone_div` | 17 | Core division and multiplication |
| `standalone_constants` | 9 | π, e, ln2, Euler-Mascheroni, sqrt constants |
| `standalone_exp_log` | 54 | exp, log, expm1, log1p, pow, log2, log10 |
| `standalone_trig` | 90 | sin, cos, tan, atan, atan2, asin, acos |
| `standalone_hyperbolic` | 57 | sinh, cosh, tanh |
| `run_validation` | 109 | Fortran reference data: add, sub, mul, div, sqrt |

**Total: 401 tests, 0 failures**

## Notes

- Build and test scripts are located in `mpfun_testarea/` (parent of the repo)
- The Kokkos install is named to reflect GPU type: `kokkos_a100/`; the build directory is `build_a100_cuda/` — this convention allows multiple GPU builds to coexist
- Cobalt output files go to the shared home directory (not `/tmp`, which is node-local)
- The `gpu_a100` queue provides a single 40 GB A100-PCIE; walltime limit is sufficient at 60 min for builds and 30 min for tests
- nvcc (CUDA 12.3) does not support C++20 generic lambdas (`[&]<int N>()`); use C++17 throughout
