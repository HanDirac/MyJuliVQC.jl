# Gradient-Descent Example

This example shows how to integrate a variational quantum circuit built with
**MyJuliVQC** into a simple **gradient-descent** training loop.  
We optimize the parameters of a VQC to **minimize the expectation value of a
Hamiltonian**, which is the standard setting in VQE-type algorithms.

Unlike the original JuliVQC example that uses `Flux.Optimise`, here we use a
plain hand-written gradient descent loop driven by `MyJuliVQC.gradient`.  
This makes the example self-contained and easy to run without extra dependencies.

```julia
using MyJuliVQC

const gradientMJVQC = MyJuliVQC.gradient

# -----------------------------
# 1. Problem setup
# -----------------------------
L      = 3                    # number of qubits
depth  = 2                    # circuit depth
ψ0     = StateVector(L)       # initial state |000⟩

# Simple 1D Heisenberg Hamiltonian as QubitsOperator:
#   H = ∑_i hz * Z_i + ∑_i J * (X_i X_{i+1} + Y_i Y_{i+1} + Z_i Z_{i+1})
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

H = heisenberg_1d(L)     # user-defined Hamiltonian as QubitsOperator

# -----------------------------
# 2. Build a variational circuit
# -----------------------------
circuit = QCircuit()

# First layer: local rotations
for i in 1:L
    push!(circuit, RzGate(i, rand(); isparas = true))
    push!(circuit, RyGate(i, rand(); isparas = true))
    push!(circuit, RzGate(i, rand(); isparas = true))
end

# Repeated entangling + rotations
for l in 1:depth
    # entangling chain
    for i in 1:L-1
        push!(circuit, CNOTGate(i, i+1))
    end
    # another layer of local rotations
    for i in 1:L
        push!(circuit, RzGate(i, rand(); isparas = true))
        push!(circuit, RxGate(i, rand(); isparas = true))
        push!(circuit, RzGate(i, rand(); isparas = true))
    end
end

# -----------------------------
# 3. Wrap loss as expectation value
#    loss(circ) = Re⟨ψ0| C(θ)† H C(θ) |ψ0⟩
# -----------------------------
loss_obj = LossExpectationRealSV(H, ψ0)

# Extract initial parameter vector
θ = active_parameters(circuit)

# -----------------------------
# 4. Simple gradient-descent loop
# -----------------------------
η      = 0.01    # learning rate
epochs = 10      # number of gradient steps

for epoch in 1:epochs
    # gradient of loss w.r.t. all active parameters in `circuit`
    gθ = gradientMJVQC(loss_obj, circuit)

    # gradient-descent update (in-place on θ)
    θ .-= η .* gθ

    # write updated parameters back into the circuit
    reset_parameters!(circuit, θ)

    # monitor current loss
    current_loss = loss_obj(circuit)
    println("Epoch $epoch: loss = $current_loss")
end
```

# Explanation of the Workflow

## 1. Circuit Construction

We manually build a parameterized quantum circuit (`QCircuit`) using  
`RzGate`, `RyGate`, `RxGate`, and `CNOTGate`.  
Gates flagged with `isparas = true` contribute entries to  
`active_parameters(circuit)` which determines the ordering of the parameter  
vector `θ`.

## 2. Hamiltonian and Loss

The Hamiltonian `H` is encoded as a `QubitsOperator`, here constructed by  
`heisenberg_1d`.

The loss is the real part of the expectation value

`⟨ψ0 | C(θ)† H C(θ) | ψ0⟩` ,

wrapped in a `LossExpectationRealSV(H, ψ0)` object.

Evaluating `loss_obj(circuit)` computes this scalar.

## 3. Gradients via `MyJuliVQC.gradient`

`MyJuliVQC.gradient` provides specialized gradient rules for  
`LossExpectationRealSV` / `LossExpectationRealDM`, internally performing a  
hand-written backward pass combined with Zygote.

Calling

```julia
gθ = gradientMJVQC(loss_obj, circuit)
```


returns the gradient vector matching the ordering of `active_parameters(circuit)`.

## 4. Hand-Written Gradient Descent

We perform simple gradient descent:

```julia
θ .-= η .* gθ
reset_parameters!(circuit, θ)
```


which updates the parameter vector and writes it back into the circuit.

---

# Notes on Future Integration with Flux / Optimisers

- This example intentionally avoids using `Flux.jl` or `Optimisers.jl` so that  
  it runs in a clean, dependency-free manner.
- The design (parameter vectors + gradients as plain `Vector{Float64}`) already  
  makes MyJuliVQC *optimizer-friendly*.
- **Future versions of MyJuliVQC may provide deeper native integration** with  
  Flux / Optimisers (e.g. ADAM or other ML optimizers), enabling smoother  
  hybrid QML–VQC workflows.


