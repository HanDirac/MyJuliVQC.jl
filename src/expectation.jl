using LinearAlgebra
using ..MyJuliVQC: QubitsOperator, QubitsTerm, DensityMatrix, nqubits
using ..MyJuliVQC: simplify!, positions, oplist, coeff, apply, QCircuit, QuantumGate, StateVector
using Zygote
#Zygote.@nograd flattened_ops
#Zygote.@nograd MyJuliVQC._flatten_ops!
Zygote.@nograd simplify!


# Internal utility for this file only: perform Kronecker products
# on several 2×2 local matrices in the specified order.
# Convention: `mats` corresponds one-to-one with the positions of a QubitsTerm, with length m ≥ 1.
_kron_many(mats::AbstractVector{<:AbstractMatrix}) = begin
    M = Matrix(mats[1])
    @inbounds for k in 2:length(mats)
        M = kron(mats[k], M)
    end
    return M
end

using Base.Threads
"""
    expectation(op::QubitsOperator, ψ::StateVector) -> Number

Compute ⟨ψ|op|ψ⟩ on a **pure state** `ψ`.  
The implementation follows these exact steps:
1) Make a copy: `op_copy = copy(op)`, then call `simplify!` to merge similar terms;
2) For each `term`:
   - Perform the Kronecker product of its local 2×2 matrices to construct an **m-qubit** `QuantumGate`;
     *Note*: here `validate=false` allows non-unitary gates (though Pauli strings are typically unitary);
   - Wrap the gate as `QCircuit(gate)` and call `apply(circ, ψ)` to obtain `|φ⟩ = term|ψ⟩`;
   - Compute `dot(ψ.data, φ.data)` (i.e. ⟨ψ|φ⟩), then multiply by the term coefficient `coeff(term)`;
3) Sum all contributions and return the result.

Does **not** materialize a full 2^n×2^n matrix.
"""
function expectation(op::QubitsOperator, ψ::StateVector)
    op_copy = copy(op); simplify!(op_copy)

    Tψ = eltype(ψ)
    partial = Vector{ComplexF64}(undef, nthreads())
    @inbounds for t in 1:length(partial)
        partial[t] = 0.0 + 0.0im
    end 

    @threads for ti in 1:length(op_copy)
        term = op_copy[ti]
        locs = Matrix{Tψ}.(oplist(term))
        Uloc = _kron_many(locs)
        φ = apply_ad(QCircuit([QuantumGate(positions(term), Uloc; validate=false)]), ψ)
        # Thread-local accumulator
        partial[threadid()] += coeff(term) * dot(ψ.data, φ.data)
    end

    s = 0.0 + 0.0im
    @inbounds for v in partial
        s += v
    end
    return s
end

# Direct trace extraction from vec(ρ):
# Tr(A) = sum(diag(A)) corresponds to indices 1:(d+1):d^2 in column-major order.
_trace_from_vec(v::AbstractVector, n::Int) = begin
    d = 1 << n
    @inbounds begin
        s = zero(eltype(v))
        idx = 1
        step = d + 1
        last = d*d
        while idx <= last
            s += v[idx]
            idx += step
        end
        return s
    end
end

# Operator application utility: vec(O ρ) = (I ⊗ O) vec(ρ)
function _apply_operator!(data::Vector{Tρ}, n::Int, pos::Vector{Int}, O::AbstractMatrix{S}) where {Tρ,S}
    pos2 = vcat(pos, n .+ pos)
    Oket = Matrix{Tρ}(O)
    Obra = Matrix{Tρ}(I, size(O,1), size(O,2))
    O2   = kron(Obra, Oket)
    data .= _apply_gate_ad(data, pos2, O2)  # modify data in place
    return nothing
end

function _apply_operator_ad(data::Vector{Tρ}, n::Int, pos::Vector{Int}, O::AbstractMatrix{S}) where {Tρ,S}
    # Similar to _apply_operator!, but does not modify data in place — returns a new array instead
    pos2 = vcat(pos, n .+ pos)
    Oket = Matrix{Tρ}(O)
    Obra = Matrix{Tρ}(I, size(O,1), size(O,2))
    O2   = kron(Obra, Oket)
    # Instead of using _apply_gate!(...) (in-place), call the small-block algorithm from apply_ad
    return _apply_gate_ad(data, pos2, O2)  # return new vector
end


"""
    expectation(op::QubitsOperator, ρ::DensityMatrix) -> Number

Compute Tr(op * ρ) on a **mixed state** `ρ`.  
Implementation steps:
1) `op_copy = copy(op)`; call `simplify!` to merge similar terms;
2) 对每个 `term`：
   - Construct the local operator `O = ⊗ᵢ Oᵢ` (aligned one-to-one with `positions(term)`, lowest qubit on the right);
   - Apply `(I ⊗ O)` to `vec(ρ)` to obtain `vec(O ρ)` (via `_apply_operator_ad`);
   - Directly compute the trace in vectorized form: `Tr(O ρ)` equals the sum of diagonal elements (stride `d+1`);
   - Multiply by `coeff(term)` and accumulate the result.
"""
function expectation(op::QubitsOperator, ρ::DensityMatrix)
    op_copy = copy(op); simplify!(op_copy)

    n = nqubits(ρ)
    
    Tψ = eltype(ρ)
    partial = Vector{ComplexF64}(undef, nthreads())
    @inbounds for t in 1:length(partial)
        partial[t] = 0.0 + 0.0im
    end

    @threads for ti in 1:length(op_copy)
        term = op_copy[ti]
        # 1) Local Kronecker product (align element type with that of ρ)
        locs = Matrix{Tψ}.(oplist(term))
        Oloc = _kron_many(locs)             # 2^m × 2^m

        # 2) Apply (I ⊗ Oloc) to vec(ρ)
        tmp = _apply_operator_ad(ρ.data, n, positions(term), Oloc)                  

        # 3) Compute Tr(O ρ) (sum diagonals in vectorized form)
        tr_val = _trace_from_vec(tmp, n)

        # 4) Multiply by coefficient and accumulate
        # Thread-local accumulator
        partial[threadid()] += coeff(term) * _trace_from_vec(tmp, n)
    end
    s = 0.0 + 0.0im
    @inbounds for v in partial
        s += v
    end
    return s
end