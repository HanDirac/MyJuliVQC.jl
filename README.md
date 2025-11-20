# MyJuliVQC

[![Build Status](https://github.com/HanDirac/MyJuliVQC.jl/actions/workflows/CI.yml/badge.svg?branch=master)](https://github.com/HanDirac/MyJuliVQC.jl/actions/workflows/CI.yml?query=branch%3Amaster)
[![Coverage](https://codecov.io/gh/HanDirac/MyJuliVQC.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/HanDirac/MyJuliVQC.jl)

## Threading Control (Added in This Fork)

⚠️ **Note:** This feature is **not part of the original MyJuliVQC** package.  
It was added in this fork (`HanDirac/MyJuliVQC`) to give users more control over threading behavior.

By default, `MyJuliVQC` will try to parallelize gate application (`apply!` / `apply_ad`) across independent slices of the state vector.  
On some machines this may conflict with BLAS internal threading, or introduce unnecessary overhead for small circuits.

You can configure the threading policy using:

```julia
using MyJuliVQC

# Example: disable outer threading, lower threshold, fix BLAS threads
MyJuliVQC.set_threading!(;
    outer_threads=false,      # whether to use outer @threads parallelism (default: true)
    dloc_threshold=8,         # gates with dloc <= threshold may trigger outer threading
    blas_threads=1            # optionally set BLAS threads (if supported on your system)
)

## Numeric Type Support (Float32, Float64, ComplexF32, ComplexF64) 

# Example

using MyJuliVQC, LinearAlgebra

# 4-qubit random normalized state (single precision)
ψ32 = StateVector(rand(ComplexF32, 16)); ψ32 ./= norm(ψ32)

# Define a simple 2-qubit gate (Float32 precision)
U = [1 0 0 0;
     0 0 1 0;
     0 1 0 0;
     0 0 0 1] .|> Float32

gate = QuantumGate([1, 2], U)
circ = QCircuit(gate)

# Apply the gate
ψ′ = apply(circ, ψ32)
println(typeof(ψ′.data))  # → Vector{ComplexF32}

# Works equivalently for Float64 / ComplexF64
ψ64 = StateVector(rand(ComplexF64, 16)); ψ64 ./= norm(ψ64)
ψ′64 = apply(circ, ψ64)
println(typeof(ψ′64.data))  # → Vector{ComplexF64}

⚙️ Notes

Internal matrix multiplications automatically use the correct element type.

Mixing precisions (e.g., ComplexF32 state with Float64 gate) will promote to the higher precision as per Julia’s rules.

When benchmarking on large systems, Float32/ComplexF32 can reduce memory footprint by 50% and often improve cache performance.

⚙️ Recommended Configuration (for HPC / Multicore)
Environment Variable	Recommended Value	Description
JULIA_NUM_THREADS	32 or 64	Julia-level threading
OMP_NUM_THREADS	8–16	BLAS / MKL threading
MKL_NUM_THREADS	8–16	Optional, for Intel MKL builds
BLAS.set_num_threads(n)	Call in Julia if needed	Manual BLAS thread control

Run Julia with -J MyJuliVQC_sysimage.so to enable precompiled sysimage.

When benchmarking or deploying on HPC, use Float32 or ComplexF32 for large-scale jobs to reduce memory and improve throughput.

The helper function set_threading!(outer_threads, dloc_threshold, blas_threads)
was added by Han Hao (2025), not part of the original JuliVQC.jl package.

## License

MyJuliVQC.jl is distributed under the GPL-3.0 license.

This project contains code adapted from the JuliVQC.jl package,
which is also licensed under GPL-3.0. Therefore, this project must
also be distributed under GPL-3.0 in compliance with the license terms.
