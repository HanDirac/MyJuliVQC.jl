# Code adapted from QuantumCircuits.jl by Weiyou Liao (original paper's reference implementation)

# src/quantum_operation.jl
# Unified abstract hierarchy + default interfaces.

abstract type QuantumOperation end
abstract type QuantumPrimitive   <: QuantumOperation end   # gates, channels, measures...
abstract type QuantumHamiltonian <: QuantumOperation end   # QubitsTerm, QubitsOperator...

# ---- Default interfaces (throwing stubs to avoid silent misuse) ----

# qubit positions involved (1-based). Containers can choose to throw.
positions(op::QuantumOperation) = throw(MethodError(positions, (op,)))
ordered_positions(op::QuantumOperation) = throw(MethodError(ordered_positions, (op,)))

# small matrix for primitives; containers or Hamiltonians may throw by default
mat(op::QuantumOperation) = throw(MethodError(mat, (op,)))
ordered_mat(op::QuantumOperation) = throw(MethodError(ordered_mat, (op,)))

# total qubit count that the object _acts on_ (if applicable)
nqubits(op::QuantumOperation) = throw(MethodError(nqubits, (op,)))

# parameter plumbing (safe defaults)
parameters(op::QuantumOperation) = ()
nparameters(op::QuantumOperation) = length(parameters(op))
