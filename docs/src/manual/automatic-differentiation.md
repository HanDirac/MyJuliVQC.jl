# Automatic Differentiation
MyJuliVQC provides transparent support for automatic differentiation (AD) of
variational quantum circuits. The core idea is:

1. Construct a parameterized circuit `circ::QCircuit`.
2. Define a **loss** based on an observable expectation value.
3. Use the specialized `gradient` API to obtain derivatives with respect to
   all active parameters in the circuit.

Internally, MyJuliVQC implements a custom, hand-written backward pass for
expectation values, integrated with Zygote, for both:

- **state-vector simulations** (`StateVector`), and
- **density-matrix simulations** (`DensityMatrix`, i.e. noisy circuits).

This is built around two lightweight loss types:

- `LossExpectationRealSV` – for pure states (state vectors).
- `LossExpectationRealDM` – for mixed states (density matrices).

These types are backed by a two-state reverse-mode algorithm that avoids
materializing full `2ⁿ × 2ⁿ` matrices.

---

## High-Level Usage

### Step 1 – Build a circuit and an operator

You first construct:

- a parameterized circuit: `circ::QCircuit`, and
- a qubit operator: `op::QubitsOperator` (e.g. a Hamiltonian).

You also choose an initial state: either a `StateVector` or `DensityMatrix`.

```julia
using MyJuliVQC

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

L = 3
ψ0 = StateVector(L)       # |000⟩
H  = heisenberg_1d(L)     # user-defined helper that builds a QubitsOperator

circ = QCircuit()
for i in 1:2
    push!(circ, CNOTGate(i, i+1))
end
for i in 1:3
    push!(circ, RyGate(i, randn(); isparas=true))
end

```

### Step 2 – Wrap the loss using LossExpectationRealSV / LossExpectationRealDM

For **pure-state** simulations:

```julia
loss_obj = LossExpectationRealSV(H, ψ0)

# This object is callable: loss_obj(circ) = real(⟨ψ0| C(θ)† H C(θ) |ψ0⟩)
val = loss_obj(circ)
println("Loss value = ", val)
```

For **density-matrix** simulations (e.g. noisy circuits):

```julia
ρ0 = DensityMatrix(L)
loss_dm = LossExpectationRealDM(H, ρ0)

val_dm = loss_dm(circ)
println("Loss (density-matrix) = ", val_dm)
```

### Step 3 – Compute gradients w.r.t. circuit parameters

MyJuliVQC provides specialized overloads of `gradient`:

```julia
using MyJuliVQC
const gradientMJVQC = MyJuliVQC.gradient

grads = gradientMJVQC(loss_obj, circ)
println("Gradient vector = ", grads)
```

The result is a `Vector{Float64}` ordered consistently with:

```julia
θ = active_parameters(circ)
```

You can then perform a parameter update, for example via gradient descent:

```julia
α = 0.01
θ  = active_parameters(circ)
θ_new = θ .- α .* grads
reset_parameters!(circ, θ_new)
```

---

## Example: Variational Energy Minimization (Pure State)

Below is a more complete example of a simple variational ansatz with MyJuliVQC’s
automatic differentiation.

```julia
using MyJuliVQC

const gradientMJVQC = MyJuliVQC.gradient

# 1. problem setup
L  = 3
ψ0 = StateVector(L)
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

H  = heisenberg_1d(L; hz = 1.0, J = 1.0)

# 2. build a parameterized circuit
circ = QCircuit()
n_layers = 2

for layer in 1:n_layers
    # entangling pattern
    for i in 1:L-1
        push!(circ, CNOTGate(i, i+1))
    end
    # local rotations
    for i in 1:L
        push!(circ, RyGate(i, randn(); isparas = true))
        push!(circ, RxGate(i, randn(); isparas = true))
    end
end

# 3. wrap loss using LossExpectationRealSV
loss_obj = LossExpectationRealSV(H, ψ0)

# 4. evaluate loss and gradient
E = loss_obj(circ)
println("Initial energy = ", E)

g = gradientMJVQC(loss_obj, circ)
println("Gradient vector length = ", length(g))

# 5. gradient descent update
α = 0.05
θ  = active_parameters(circ)
θ_new = θ .- α .* g
reset_parameters!(circ, θ_new)

E_new = loss_obj(circ)
println("Energy after one GD step = ", E_new)
```

---

## Example: Variational Optimization with Noise (Density Matrix)

For noisy circuits, you only need to replace the state with a density matrix
and use `LossExpectationRealDM`:

```julia
using MyJuliVQC

const gradientMJVQC = MyJuliVQC.gradient

L  = 2
ρ0 = DensityMatrix(L)
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

circ = QCircuit()
push!(circ, HGate(1))
push!(circ, Depolarizing(1;p=0.05))
push!(circ, RxGate(1, 0.3; isparas = true))
push!(circ, CNOTGate(1, 2))

loss_dm = LossExpectationRealDM(H, ρ0)

E = loss_dm(circ)
println("Noisy loss = ", E)

g = gradientMJVQC(loss_dm, circ)
println("Gradient (noisy case) = ", g)
```
---

## Practical Notes and Limitations

- **Currently, MyJuliVQC only supports gradients for losses constructed using**
  `LossExpectationRealSV` **and** `LossExpectationRealDM`.

- These two loss types correspond to expectation-based objectives of the form  
  `Re(⟨ψ| C(θ)† H C(θ) |ψ⟩)` for pure states and  
  `Re(Tr[ H C(θ) ρ C(θ)† ])` for density matrices.

- Any *other* form of loss function (e.g., fidelity loss, distance-based loss,
  norm penalties, custom metric functions, etc.) is **not yet supported** by  
  `MyJuliVQC.gradient`. Attempting to call `gradient(loss, circ)` on such losses
  will not work.

- Full AD support for arbitrary loss functions may be added in a future version
of MyJuliVQC.

---

## Summary

- MyJuliVQC provides custom AD support for expectation values using
  `LossExpectationRealSV` and `LossExpectationRealDM`.
- The public API is:

  ```julia
  loss_obj = LossExpectationRealSV(op, ψ0)
  grads    = gradient(loss_obj, circ)
  ```

  (and analogously for `LossExpectationRealDM`).

- Gradients are returned as a `Vector{Float64}` aligned with
  `active_parameters(circ)` and can be directly used in gradient-based
  optimization loops for VQE, variational classifiers, and other hybrid
  quantum–classical algorithms.
