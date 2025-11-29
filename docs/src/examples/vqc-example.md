# VQC-example

A variational quantum circuit (VQC) consists of quantum gates with parameters that can be optimized. The code below demonstrates how to build one step by step.

```julia
using MyJuliVQC

const gradientMJVQC = MyJuliVQC.gradient

# -----------------------------
# 1. Circuit parameters
# -----------------------------
L     = 3          # number of qubits
depth = 2          # circuit depth

# -----------------------------
# 2. Build a variational circuit
# -----------------------------
circuit = QCircuit()

# First layer: local rotations
for i in 1:L
    push!(circuit, RzGate(i, rand();  isparas = true))
    push!(circuit, RyGate(i, rand();  isparas = true))
    push!(circuit, RzGate(i, rand();  isparas = true))
end

# Repeated entangling + rotations
for l in 1:depth
    # entangling CNOT chain
    for i in 1:L-1
        push!(circuit, CNOTGate(i, i+1))
    end
    # another layer of local rotations
    for i in 1:L
        push!(circuit, RzGate(i, rand();  isparas = true))
        push!(circuit, RxGate(i, rand();  isparas = true))
        push!(circuit, RzGate(i, rand();  isparas = true))
    end
end

# -----------------------------
# 3. Initial state & Hamiltonian
# -----------------------------
ψ0 = StateVector(L)              # |000⟩
function heisenberg_1d(L; hz = 1.0, J = 1.0)
    terms = QubitsTerm[]

    # Local Z fields
    for i in 1:L
        push!(terms, QubitsTerm(i => "Z"; coeff = hz))
    end

    # Nearest-neighbor interactions: X X + Y Y + Z Z
    for i in 1:L-1
        push!(terms, QubitsTerm(i => "X", i+1 => "X"; coeff = J))
        push!(terms, QubitsTerm(i => "Y", i+1 => "Y"; coeff = J))
        push!(terms, QubitsTerm(i => "Z", i+1 => "Z"; coeff = J))
    end

    return QubitsOperator(terms)
end
H  = heisenberg_1d(L)            

# -----------------------------
# 4. Wrap loss as an expectation value
#    loss(circ) = Re⟨ψ0| C(θ)† H C(θ) |ψ0⟩
# -----------------------------
loss_obj = LossExpectationRealSV(H, ψ0)
E = loss_obj(circuit)
println("Initial loss (energy) = ", E)

# -----------------------------
# 5. Extract parameters, compute gradient, update
# -----------------------------
θ  = active_parameters(circuit)
println("Number of parameters = ", length(θ))

grads = gradientMJVQC(loss_obj, circuit)
println("Gradient vector length = ", length(grads))

# Gradient descent
α      = 0.01
θ_new  = θ .- α .* grads
reset_parameters!(circuit, θ_new)

E_new = loss_obj(circuit)
println("Loss after one gradient step = ", E_new)
```
