# MyJuliVQC

[![Build Status](https://github.com/HanDirac/MyJuliVQC.jl/actions/workflows/CI.yml/badge.svg?branch=master)](https://github.com/HanDirac/MyJuliVQC.jl/actions/workflows/CI.yml?query=branch%3Amaster)
[![Coverage](https://codecov.io/gh/HanDirac/MyJuliVQC.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/HanDirac/MyJuliVQC.jl)

**MyJuliVQC.jl** is an independent Julia implementation of variational quantum circuit (VQC) simulation methods inspired by the original research work and design philosophy of JuliVQC by Weiyou Liao.  
This package reimplements key algorithms based on the corresponding scientific publication:
**[JuliVQC: an Efficient Variational Quantum Circuit Simulator for Near-Term Quantum Algorithms](https://link.springer.com/article/10.1140/epjs/s11734-025-01628-6))**
while adopting a clean and modular architecture tailored for clarity, extensibility, and educational use.

Although this project is technically independent from **JuliVQC.jl**, it aims to reproduce a similar set of functionalities in a self-contained way, making it suitable for:

- building and simulating parameterized quantum circuits,
- exploring variational quantum algorithms,
- prototyping quantum machine learning models,
- and understanding how circuit-level simulation works under the hood in Julia.

> Note: This repository includes several functions adapted from JuliVQC's open-source implementation when the papers did not explicitly specify algorithmic details.  
> Following GPL-3.0 licensing requirements, this project is also released under GPL-3.0.

---

## ✨ Features

MyJuliVQC provides a clean, modular, and research-oriented implementation of variational quantum circuit simulation.  
Although designed for clarity and educational value, the package covers a comprehensive set of essential tools for near-term quantum algorithm prototyping.

### 🚀 Core Features

- **State-Vector Simulation Engine**  
  Fully supports quantum state initialization, evolution, and measurement with user-defined or built-in quantum gates.

- **Modular Quantum Circuit Construction**  
  Flexible `QCircuit` structure with intuitive APIs for assembling parameterized quantum circuits.

- **Support for Parameterized Gates (VQC-ready)**  
  Native support for variational parameters, enabling VQE, QML models, and other hybrid quantum-classical algorithms.

- **Hamiltonian & Observable Utilities**  
  Built-in tools for constructing qubit Hamiltonians, evaluating expectation values, and computing gradient-based losses.

- **Readable, Transparent Codebase**  
  MyJuliVQC focuses on correctness and clarity, making it ideal for learning, research demonstrations, and rapid prototyping.

- **Algorithm-Guided Implementation**  
  Core algorithms follow the descriptions in variational quantum simulation literature, ensuring theoretical fidelity.

- **Self-Contained Execution**  
  No external simulator required — MyJuliVQC runs using a standalone Julia implementation without depending on external quantum SDKs.

### ⚡ Performance Notice

While MyJuliVQC prioritizes clarity and modularity over raw performance, its design remains fully compatible with VQC workflows.  
Even with modest speed, the package provides a reliable foundation for implementing, analyzing, and experimenting with variational algorithms.



## 📦 Installation

Install MyJuliVQC directly from GitHub:

```julia
julia> ] add https://github.com/HanDirac/MyJuliVQC.jl
julia> using MyJuliVQC
```

---

## 🚀 Quick Examples

Here are a few minimal examples demonstrating how to use **MyJuliVQC** to construct and simulate variational quantum circuits.  
The API is designed to be intuitive and lightweight, making it easy to explore both basic quantum operations and variational workflows.

---

### Example 1: Preparing a Bell State

The following script initializes a 2-qubit state, applies a simple quantum circuit, and performs a measurement:

```julia
using MyJuliVQC

# Create a 2-qubit state |00⟩
ψ = StateVector(2)

# Build a quantum circuit
circuit = QCircuit()
push!(circuit, HGate(1))          # Hadamard on qubit 1
push!(circuit, CNOTGate(1, 2))    # Controlled-NOT from qubit 1 to 2

# Apply the circuit to the state
apply!(circuit, ψ)

# Inspect the resulting state
println("Final state vector:")
println(ψ)

# Measure qubit 1
outcome, prob = measure!(ψ, 1)
println("Measured qubit 1 → $outcome   (probability = $prob)")
```

---

### Example 2: A Simple Variational Circuit

This example shows how to construct a parameterized circuit and evaluate a loss function:

```julia
using MyJuliVQC

L = 3                        # number of qubits
ψ0 = StateVector(L)          # initial state

# Build a simple variational circuit
circuit = QCircuit()

# Layer 1: single-qubit rotations
for i in 1:L
    push!(circuit, RzGate(i, rand(), isparas=true))
    push!(circuit, RyGate(i, rand(), isparas=true))
    push!(circuit, RzGate(i, rand(), isparas=true))
end

# Layer 2: entangling + rotations
for i in 1:L-1
    push!(circuit, CNOTGate(i, i+1))
end
for i in 1:L
    push!(circuit, RxGate(i, rand(), isparas=true))
end

# Define a target state (example purposes only)
target = StateVector(rand(ComplexF64, 2^L))
normalize!(target)

# Define a simple loss function
loss(c) = distance(target, c * ψ0)

println("Loss value: ", loss(circuit))
```

---

You can use these examples as a starting point for building VQE solvers, variational classifiers, and other hybrid quantum-classical applications.

## 📚 Documentation

Documentation is built using Documenter.jl.

You can view the online documentation at:

https://<your-username>.github.io/MyJuliVQC.jl

(To enable this, make sure GitHub Pages and Documenter CI workflow are set up.)

---

## 🧪 Testing

Tests are located in:

```
test/
```

Run using:

```julia
julia> ] test MyJuliVQC
```

---

## 📂 Project Structure

```
MyJuliVQC.jl/
  src/
  test/
  docs/
  examples/   (optional)
  Project.toml
  LICENSE
  README.md
```

---

## 🔬 Scientific Background

MyJuliVQC is developed based on the variational quantum simulation framework described in the following works:

- (Add your target papers here)
- (Include papers corresponding to JuliVQC if needed)

This package reimplements the algorithms described therein with additional engineering choices for clarity and modularity.

---

## 👤 Authors

**Han Hao (韩昊)**  
Tang Aoqing Honor Program in Science  
Jilin University  
Email: 515673679@qq.com

---

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
