# MPFUN2020 → Kokkos Port: Project Plan

## Executive Summary

Port David Bailey's MPFUN2020 arbitrary-precision floating-point library from Fortran-90 to a C++/Kokkos implementation that runs on CPUs, NVIDIA GPUs (CUDA), AMD GPUs (HIP), and Intel GPUs (SYCL).

**Source:** MPFUN20-Fort v33 (~39k lines Fortran-90)  
**Target:** Header-only or minimal-library C++17 with Kokkos parallel constructs  
**Goal:** `KOKKOS_INLINE_FUNCTION`-compatible arbitrary precision arithmetic

---

## Source Code Analysis

### Module Structure (fortran-var1)

| File | Lines | Purpose |
|------|-------|---------|
| `mpfuna.f90` | 986 | Global parameters, constants (log2, pi, egamma precomputed) |
| `mpfunb.f90` | 3,590 | **Core arithmetic**: add, sub, mul, div, sqrt, compare |
| `mpfunc.f90` | 963 | I/O, conversion (DP↔MPR, string↔MPR) |
| `mpfund.f90` | 1,743 | Transcendentals: exp, log, sin, cos, tan, atan, etc. |
| `mpfune.f90` | 9,412 | Special functions: gamma, zeta, Bessel, hypergeometric, etc. |
| `mpfunf.f90` | 65 | User precision parameter (`mpipl`) |
| `mpfung1.f90` | 4,922 | High-level wrappers (variant 1) |
| `mpfunh1.f90` | 4,886 | Additional high-level wrappers |
| `mpmodule.f90` | 92 | Module umbrella |

### Key Design Decisions in Original

1. **Representation**: Array of 64-bit integers, 60 bits used per word
   - `mpnbt = 60` (bits per mantissa word)
   - `mpbdx = 2^60` (radix)
   - First words are metadata: `a(0)` = allocated size, `a(1)` = working precision, `a(2)` = sign+length, `a(3)` = exponent

2. **Precision**: Configurable at compile time via `mpipl` (digits), converted to `mpwds` (words)
   - `mpwds = int(mpipl / mpdpw + 2)` where `mpdpw ≈ 18.06`

3. **FFT multiplication**: Used for very high precision (>~20k digits)

4. **Thread safety**: No global mutable state in computation paths

---

## Port Architecture

### Design Principles

1. **Fixed-size storage** for GPU compatibility (no dynamic allocation in kernels)
2. **Template-based precision** at compile time: `MPFloat<WORDS>`
3. **Kokkos-compatible**: All core functions marked `KOKKOS_INLINE_FUNCTION`
4. **Expression templates** for operator fusion (optional, phase 2)
5. **No STL containers** in device code; use Kokkos Views where needed

### Proposed C++ Structure

```
kokkos-mpfun/
├── include/
│   └── mpfun/
│       ├── mpfun.hpp           # Single-include header
│       ├── config.hpp          # Compile-time config (precision, features)
│       ├── mp_float.hpp        # Core MPFloat<N> type
│       ├── mp_complex.hpp      # MPComplex<N> type
│       ├── core/
│       │   ├── representation.hpp  # Data layout, accessors
│       │   ├── add_sub.hpp         # Addition, subtraction
│       │   ├── mul.hpp             # Multiplication (schoolbook + FFT)
│       │   ├── div.hpp             # Division
│       │   ├── sqrt.hpp            # Square root
│       │   └── compare.hpp         # Comparison operators
│       ├── transcendental/
│       │   ├── exp_log.hpp
│       │   ├── trig.hpp
│       │   └── hyperbolic.hpp
│       ├── special/              # Phase 2+
│       │   ├── gamma.hpp
│       │   ├── bessel.hpp
│       │   └── ...
│       ├── io/                   # Host-only
│       │   ├── string_conv.hpp
│       │   └── format.hpp
│       └── constants/
│           ├── pi.hpp
│           ├── log2.hpp
│           └── egamma.hpp
├── tests/
│   ├── unit/
│   └── validation/             # Compare against MPFUN reference
├── benchmarks/
└── examples/
```

---

## Implementation Phases

### Phase 1: Core Arithmetic (MVP)
**Goal**: Addition, subtraction, multiplication, division, comparison, sqrt  
**Deliverable**: Can compute `(a + b) * c / d` at arbitrary precision on GPU

| Component | Fortran Source | Est. Effort |
|-----------|---------------|-------------|
| Representation | `mpfuna.f90` (params) | 1 day |
| Add/Sub | `mpfunb.f90:mpadd,mpsub,mpnorm` | 2 days |
| Multiply | `mpfunb.f90:mpmul,mpmulx` | 3 days |
| Divide | `mpfunb.f90:mpdiv,mpdivd` | 2 days |
| Compare | `mpfunb.f90:mpcpr` | 0.5 days |
| Sqrt | `mpfunb.f90:mpsqrt` | 1 day |
| **Phase 1 Total** | | **~10 days** |

### Phase 2: Transcendentals
**Goal**: exp, log, sin, cos, tan, atan, sinh, cosh, tanh

| Component | Fortran Source | Est. Effort |
|-----------|---------------|-------------|
| Exp | `mpfund.f90:mpexp` | 1.5 days |
| Log | `mpfund.f90:mplog` | 1.5 days |
| Sin/Cos | `mpfund.f90:mpcssn` | 2 days |
| Tan | `mpfund.f90:mptan` | 0.5 days |
| Atan/Atan2 | `mpfund.f90:mpang` | 1 day |
| Hyperbolics | `mpfund.f90:mpcssh` | 1 day |
| **Phase 2 Total** | | **~8 days** |

### Phase 3: Complex Numbers
**Goal**: Complex arithmetic mirroring MPComplex

| Component | Est. Effort |
|-----------|-------------|
| Complex type | 1 day |
| Complex arithmetic | 2 days |
| Complex transcendentals | 2 days |
| **Phase 3 Total** | **~5 days** |

### Phase 4: Special Functions (Optional)
**Goal**: Gamma, Bessel, Zeta, etc.  
**Estimate**: 2-4 weeks depending on scope

### Phase 5: Optimization
- FFT-based multiplication for high precision
- Expression templates
- Memory layout tuning (AoS vs SoA)
- Architecture-specific tuning

---

## Agent Orchestration Strategy

### Recommended: 3-4 Parallel Agents

Given the modular structure and clear interfaces, we can parallelize effectively:

```
┌─────────────────────────────────────────────────────────────────┐
│                     COORDINATOR (Main Session)                   │
│  - Overall architecture decisions                                │
│  - Interface definitions                                         │
│  - Integration and conflict resolution                          │
│  - Final review and testing                                      │
└─────────────────────────────────────────────────────────────────┘
                              │
         ┌────────────────────┼────────────────────┐
         ▼                    ▼                    ▼
┌─────────────────┐  ┌─────────────────┐  ┌─────────────────┐
│   AGENT 1       │  │   AGENT 2       │  │   AGENT 3       │
│   Core Arith    │  │   Transcendental│  │   Testing/Infra │
├─────────────────┤  ├─────────────────┤  ├─────────────────┤
│ • Representation│  │ • exp/log       │  │ • Test harness  │
│ • add/sub       │  │ • sin/cos/tan   │  │ • Validation    │
│ • mul (basic)   │  │ • atan          │  │ • Benchmarks    │
│ • div           │  │ • hyperbolics   │  │ • CI setup      │
│ • sqrt          │  │                 │  │ • Examples      │
│ • compare       │  │ (waits for A1   │  │                 │
│                 │  │  core types)    │  │                 │
└─────────────────┘  └─────────────────┘  └─────────────────┘
```

### Why 3-4 Agents (not more)?

1. **Dependency chain**: Transcendentals depend on core arithmetic; can't parallelize infinitely
2. **Interface coherence**: More agents = more integration overhead
3. **Code review bandwidth**: You need to review and approve merges
4. **Diminishing returns**: The Fortran is ~40k lines, but much is comments/wrappers

### Agent Roles

#### Agent 1: Core Arithmetic Lead
- Owns `mp_float.hpp`, `representation.hpp`, `core/*.hpp`
- Defines the fundamental data type and memory layout
- Implements add, sub, mul, div, sqrt, compare
- **Critical path**: Others depend on this

#### Agent 2: Transcendentals
- Owns `transcendental/*.hpp`
- Starts with stub interface, implements once Agent 1 delivers core
- Can work on algorithm design while waiting

#### Agent 3: Test Infrastructure & Validation
- Sets up CMake/build system
- Creates test harness comparing against reference Fortran
- Writes validation tests (compute known values, compare)
- Sets up benchmarking framework

#### Agent 4 (Optional): Complex Numbers
- Can work in parallel once core real arithmetic exists
- Owns `mp_complex.hpp`

### Communication Protocol

1. **Interface contracts first**: Before coding, agents agree on:
   - `MPFloat<N>` type signature
   - Memory layout (`data_[N+4]` array structure)
   - Function signatures for core operations

2. **Shared header for types**: `representation.hpp` written first by Agent 1, reviewed by all

3. **Integration points**:
   - After Agent 1 completes basic add/mul: Agent 2 can start
   - After test harness exists: All agents run validation

### No Deep Nesting

Nested sub-agents add coordination overhead without clear benefit here. The work decomposes cleanly into 3-4 parallel streams. If an agent needs to delegate subtasks, they can do so internally without spawning persistent child agents.

---

## Technical Decisions

### 1. Template Precision vs Runtime Precision

**Decision**: Template-based (compile-time) precision

```cpp
template <int WORDS>
class MPFloat {
    int64_t data_[WORDS + 4];  // +4 for metadata
    // ...
};

using MP100 = MPFloat<6>;    // ~100 digits
using MP1000 = MPFloat<57>;  // ~1000 digits
```

**Rationale**:
- GPU-compatible (no dynamic allocation)
- Better optimization opportunities
- Matches MPFUN's compile-time `mpipl` parameter

### 2. Memory Layout

Following MPFUN's proven layout:
```
data_[0] = allocated_size
data_[1] = working_precision (words)
data_[2] = sign * mantissa_length
data_[3] = exponent
data_[4...] = mantissa words (60 bits used per 64-bit word)
```

### 3. Kokkos Integration

```cpp
template <int WORDS>
class MPFloat {
public:
    KOKKOS_INLINE_FUNCTION MPFloat() = default;
    KOKKOS_INLINE_FUNCTION MPFloat(double d);
    
    KOKKOS_INLINE_FUNCTION MPFloat operator+(const MPFloat& rhs) const;
    // ...
    
private:
    int64_t data_[WORDS + 4];
};
```

### 4. FFT Multiplication

- **Phase 1**: Schoolbook O(n²) multiplication (simpler, works for moderate precision)
- **Phase 5**: FFT-based multiplication for >1000 digits
  - Will need Kokkos-compatible FFT (cuFFT/rocFFT wrappers or custom)

---

## Validation Strategy

1. **Unit tests**: Each function against known values
2. **Cross-validation**: Run same computation in Fortran MPFUN and C++ port, compare results
3. **Convergence tests**: Series expansions that should converge to known constants
4. **Stress tests**: Large parallel workloads on GPU

### Reference Test Cases (from MPFUN)

- `testmpfun.f90`: Comprehensive function tests
- `tpslq1.f90`: PSLQ algorithm (integer relation detection)
- `tquad.f90`: Numerical integration

---

## Risk Assessment

| Risk | Likelihood | Impact | Mitigation |
|------|------------|--------|------------|
| FFT portability | Medium | High | Start with schoolbook mul; FFT in phase 5 |
| Precision edge cases | Medium | Medium | Extensive validation against Fortran reference |
| Performance regression | Medium | Medium | Benchmark early, profile GPU kernels |
| Kokkos version issues | Low | Medium | Target Kokkos 4.x, document requirements |

---

## Timeline Estimate

| Phase | Duration | Parallelizable? |
|-------|----------|-----------------|
| Phase 1: Core | 2 weeks | Partially (A1 + A3) |
| Phase 2: Transcendentals | 1.5 weeks | Yes (A2 parallel to late A1) |
| Phase 3: Complex | 1 week | Yes (parallel to Phase 2) |
| Integration/Polish | 1 week | No |
| **Total** | **5-6 weeks** | |

With 3 parallel agents working effectively, could compress to **3-4 weeks** for phases 1-3.

---

## Next Steps

1. **Define interface contracts** (`representation.hpp`, `MPFloat` API)
2. **Set up repository structure** and build system
3. **Spawn agents** with clear task assignments
4. **Begin Phase 1** implementation

---

## Appendix: Key Fortran Routines to Port

### Critical Path (Phase 1)
- `mpadd`, `mpsub` - Addition/subtraction
- `mpmul`, `mpmulx` - Multiplication
- `mpdiv`, `mpdivd` - Division
- `mpsqrt` - Square root
- `mpcpr` - Comparison
- `mpnorm` - Normalization
- `mpeq` - Assignment/copy

### Phase 2 Dependencies
- `mpexp` - Exponential (uses mul, add)
- `mplog` - Logarithm (uses div, mul, add)
- `mpcssn` - Sin/Cos (uses mul, add, sub)
- `mpang` - Arctangent (uses div, mul)
- `mpcssh` - Sinh/Cosh (uses exp)
