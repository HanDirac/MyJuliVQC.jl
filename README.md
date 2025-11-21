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
const gradientMJVQC = MyJuliVQC.gradient

# 1. initial state
L  = 3                      # number of qubits
ψ0 = StateVector(L)         # |000⟩

# 2. construct a simple variational circuit
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

# 3. define a simple operator H
H = QubitsOperator([QubitsTerm(1=>"Z", 2=>"Z", 3=>"Z"; coeff=1.0)])

# 4. construct the loss
loss_obj = LossExpectationRealSV(H, ψ0)

# 5. compute the loss(expectation) with the current circuit
println("Expectation value (loss) = ", loss_obj(circuit))

# 6. use gradient to get the gradient
grads = gradientMJVQC(loss_obj, circuit)
println("Gradient from MyJuliVQC.gradient:")
println(grads)

```

---

You can use these examples as a starting point for building VQE solvers, variational classifiers, and other hybrid quantum-classical applications.

## 📚 Documentation

The documentation for **MyJuliVQC** will be available here:

👉 https://HanDirac.github.io/MyJuliVQC.jl/

---

### Local Deployment 🖥️

If you would like to build and view the documentation locally (once the docs are added),  
you may follow the typical Julia *Documenter.jl* workflow:

1. Clone this repository:
   ```bash
   git clone https://github.com/HanDirac/MyJuliVQC.jl.git
   cd MyJuliVQC.jl
   ```
   
2. Install documentation dependencies inside the docs environment:
    ```julia
     cd docs
     julia --project=.
     using Pkg
     Pkg.instantiate()
     ```

3. Build the documentation:
    ```julia
     julia make.jl
     ```

4. View the documentation:
- The built HTML file will be located at:
  ```
  docs/build/index.html
  ```
- Open this file in your browser to preview the documentation locally.

---

## 📂 Project Structure

```
MyJuliVQC.jl/
  src/
  test/
  docs/
  Project.toml
  LICENSE
  README.md
```

---


## 👤 Author

**Han Hao (郝瀚)**  
School of Physics, Jilin University  
(Part of the development was completed during my research assistantship at the Hefei National Laboratory.)  
Email: 515673679@qq.com


## Advanced Topics

Detailed documentation on threading control, numeric precision, and HPC tuning is available in the online documentation:

👉 https://HanDirac.github.io/MyJuliVQC.jl/
