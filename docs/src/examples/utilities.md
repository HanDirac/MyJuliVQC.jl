# Utilities

The following utility functions are provided by **MyJuliVQC.jl** to make it easier
to inspect and manipulate variational circuits:

```julia
active_parameters(circ::QCircuit)
reset_parameters!(circ::QCircuit, coeffs::AbstractVector{<:Number})
fuse_gates(circ::QCircuit) :: QCircuit
```

- `active_parameters(circ::QCircuit)`  
  Returns a `Vector{Float64}` containing all *active* (trainable) parameters in the
  circuit.  
  Active parameters correspond to gates constructed with `isparas = true`.  
  The order of parameters in this vector matches the internal left-to-right
  traversal of `circ`.

- `reset_parameters!(circ::QCircuit, coeffs::AbstractVector{<:Number})`  
  Overwrites all active parameters in `circ` using the entries of `coeffs`.  
  Useful for:
  - inserting a new parameter vector (e.g., produced by an optimizer),
  - scanning over parameter values,
  - evaluating the circuit at different parameter sets.  
  The length of `coeffs` must match `length(active_parameters(circ))`.

- `fuse_gates(circ::QCircuit) :: QCircuit`  
  Performs a **preliminary noiseless simplification** of the circuit.  
  If a **single-qubit gate** is adjacent to a **two-qubit gate** and acts on one
  of its qubits, it is *absorbed* into the two-qubit gate by matrix multiplication.  
  Key points:
  - Only applies to **noiseless** circuits (no `QuantumMap` allowed).
  - Absorption is **local** and **adjacency-based** (left and/or right neighbor).
  - A fused two-qubit gate becomes a general `QuantumGate` with an updated 4×4 matrix.
  - Intended as a light-weight preprocessing step before simulation or training.

These utilities support common workflows such as:
- extracting the optimization parameter vector,
- pushing updated parameters back into the circuit,
- optionally simplifying circuits prior to simulation.

They collectively make MyJuliVQC more convenient for variational algorithms,
VQE/VQC experiments, and custom optimization routines.
