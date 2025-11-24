# Quick Start

This page provides a short hands-on introduction to **MyJuliVQC**, showing how to construct, apply, and differentiate variational quantum circuits.

---
## Example 1: Preparing a Bell State

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

## Example 2: A Simple Variational Circuit

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

