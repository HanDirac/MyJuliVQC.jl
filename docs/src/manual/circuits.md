# Circuits

In **MyJuliVQC**, a quantum circuit is represented by a lightweight wrapper type `QCircuit`,
which internally stores a sequence of quantum operations.

Each element of a `QCircuit` can be:

- a (possibly parametric) **unitary gate** (e.g. `HGate`, `RxGate`, `CNOTGate`, …),
- a **quantum channel** (e.g. `AmplitudeDamping`, `Depolarizing`, or a general `QuantumMap`),
- or another `QCircuit` (allowing circuit composition via concatenation of subcircuits).

This design keeps the circuit object simple and transparent while still being flexible enough to
represent both noiseless and noisy variational circuits.

---

## Applying a Circuit to a State

After constructing a circuit, you can apply it to a quantum state using:

- `apply!(circuit, state)` – **in-place** evolution (modifies `state` directly)
- `apply(circuit, state)` – **out-of-place** evolution (returns a new state)
- `circuit * state` – equivalent to `apply(circuit, state)`

Here, `state` can be either:

- `StateVector` – representing a pure state, or
- `DensityMatrix` – representing a mixed state (needed when quantum channels are present).

For variational quantum algorithms, the out-of-place form is often convenient, since it allows you
to treat the circuit as a pure transformation and reuse the initial state.

---

## Example: Building and Running a Simple Circuit

```julia
using MyJuliVQC

# Initialize a 2-qubit |00⟩ state
state = StateVector(2)

# Construct a circuit:
#   1. H on qubit 1
#   2. Ry(π/4) on qubit 1 (non-parametric here)
#   3. CNOT with control = 1, target = 2
circuit = QCircuit([
    HGate(1),
    RyGate(1, π/4; isparas = false),
    CNOTGate(1, 2),
])

# Apply the circuit in-place
apply!(circuit, state)

# Measure qubit 2
outcome, prob = measure!(state, 2)
println("Measurement result on qubit 2 = $outcome (prob = $prob)")
```

You can also build the same circuit incrementally:

```julia
using MyJuliVQC

state = StateVector(2)
circuit = QCircuit()

push!(circuit, HGate(1))
push!(circuit, RyGate(1, π/4; isparas = false))
push!(circuit, CNOTGate(1, 2))

ψ_out = apply(circuit, state)   # out-of-place application
```

In more advanced workflows, parametric gates (with `isparas = true`) can be used together with
the gradient engine described in the **Automatic differentiation** section to implement VQE,
quantum classifiers, and other variational algorithms.