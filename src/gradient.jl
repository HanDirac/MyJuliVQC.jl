using ChainRulesCore: NoTangent
using LinearAlgebra
using Zygote
using ..MyJuliVQC:
    # circuit & primitives
    QCircuit, QuantumGate, apply, nqubits,
    qubits, matrix, active_parameters, reset_parameters!
    # states & expectation
    StateVector, DensityMatrix, QubitsOperator, expectation,
    # term helpers
    positions, oplist, coeff, kraus
using ..MyJuliVQC.Gates: ParamOp

# We need to call the internal “apply local gate on global state” core function
# (defined in apply.jl as _apply_gate!), so we reference it by fully qualified name:
const _apply_gate! = MyJuliVQC._apply_gate!

_with_params(circ, θ) = reset_parameters!(deepcopy(circ), θ)

# -----------------------------------------------------------------------------
# Utility: Apply one term of an operator H (= ∑ₜ cₜ ⊗ₖ Oₜₖ) to |ψ⟩ and return H|ψ⟩
# (without materializing a 2ⁿ×2ⁿ matrix; uses the internal _apply_gate! iteratively)
# -----------------------------------------------------------------------------

const _kron_many = MyJuliVQC._kron_many

# out = O|in>
function _apply_local_out(in_state::StateVector{Tψ},
                          pos::Vector{Int},
                          O::AbstractMatrix) where {Tψ<:Number}
    out = StateVector{Tψ}(copy(in_state.data))
    _apply_gate!(out.data, pos, O)
    return out
end

# |Ψ⟩ = H|ψ⟩
function _apply_operator_on_state(op::QubitsOperator, ψ::StateVector)
    Tψ = eltype(ψ)
    acc = zeros(Tψ, length(ψ))
    for t in op
        Uloc = _kron_many(Matrix{Tψ}.(oplist(t)))     # 2^m × 2^m
        φ = _apply_local_out(ψ, positions(t), Uloc)   # φ = (⊗ O_k) |ψ⟩
        @inbounds @simd for k in eachindex(acc)
            acc[k] += coeff(t) * φ.data[k]
        end
    end
    return StateVector{Tψ}(acc)
end

# Apply op to vec(ρ): Ψ = Ĥ * vec(ρ) = Σₜ (I ⊗ Oₜ) * vec(ρ)
const _apply_operator!_dm = MyJuliVQC._apply_operator! 

function _apply_operator_on_dmvec(op::QubitsOperator, n::Int, vρ::AbstractVector)
    T = eltype(vρ)
    acc = zeros(T, length(vρ))
    tmp = similar(acc)
    for t in op
        copyto!(tmp, vρ)
        # (I ⊗ Oloc) * tmp
        Oloc = _kron_many(Matrix{T}.(oplist(t)))
        _apply_operator!_dm(tmp, n, positions(t), Oloc)
        @inbounds @simd for k in eachindex(acc)
            acc[k] += coeff(t) * tmp[k]
        end
    end
    return acc
end

# For a parameterized operation: if it provides op.diff(params, k)::Matrix or ::QuantumGate,
# use it directly; otherwise use central difference (default eps = 1e-6).
#function _local_derivative_matrix(op, k::Int; eps::Float64=1e-6)
#    p = op.params
#    hasdiff = haskey(op, :diff) && (op.diff !== nothing)
#    if hasdiff
#        D = op.diff(p, k)
#        return D isa QuantumGate ? matrix(D) : Matrix(D)
#    else
#        p⁺ = copy(p); p⁻ = copy(p)
#        p⁺[k] += eps; p⁻[k] -= eps
#        U⁺ = op.build(p⁺); U⁻ = op.build(p⁻)
#        M⁺ = U⁺ isa QuantumGate ? matrix(U⁺) : Matrix(U⁺)
#        M⁻ = U⁻ isa QuantumGate ? matrix(U⁻) : Matrix(U⁻)
#        return (M⁺ .- M⁻) ./ (2eps)
#    end
#end

# For ParamOp (struct): no :diff field, so use central difference
function _local_derivative_matrix(op::ParamOp, k::Int; eps::Float64=1e-6)
    p = copy(op.params)
    p⁺ = copy(p); p⁻ = copy(p)
    p⁺[k] += eps
    p⁻[k] -= eps
    U⁺ = op.build(p⁺)
    U⁻ = op.build(p⁻)
    M⁺ = U⁺ isa QuantumGate ? matrix(U⁺) : Matrix(U⁺)
    M⁻ = U⁻ isa QuantumGate ? matrix(U⁻) : Matrix(U⁻)
    return (M⁺ .- M⁻) ./ (2eps)
end

# NamedTuple form (may have :diff optionally)
function _local_derivative_matrix(op::NamedTuple, k::Int; eps::Float64=1e-6)
    p = op.params
    if haskey(op, :diff) && (op.diff !== nothing)
        D = op.diff(p, k)
        return D isa QuantumGate ? matrix(D) : Matrix(D)
    else
        p⁺ = copy(p); p⁻ = copy(p)
        p⁺[k] += eps; p⁻[k] -= eps
        U⁺ = op.build(p⁺); U⁻ = op.build(p⁻)
        M⁺ = U⁺ isa QuantumGate ? matrix(U⁺) : Matrix(U⁺)
        M⁻ = U⁻ isa QuantumGate ? matrix(U⁻) : Matrix(U⁻)
        return (M⁺ .- M⁻) ./ (2eps)
    end
end

# Convert circuit (flattened) into “layers”:
# Each layer stores pos::Vector{Int}, U::Matrix, param_idx::Vector{Int}, deriv::Function(k->dU)
function _materialize_layers(c::QCircuit)
    ops = c
    layers = Vector{NamedTuple}(undef, length(ops))
    for (i, op) in enumerate(ops)
        if op isa ParamOp
            U  = op.build(op.params)
            Um = U isa QuantumGate ? matrix(U) : Matrix(U)
            pos = U isa QuantumGate ? qubits(U) : Int[1]
            idx = Int[]
            for k in eachindex(op.params)
                op.mask[k] && push!(idx, k)
            end
            layers[i] = (; pos, U=Um, param_idx=idx, deriv=(k-> _local_derivative_matrix(op,k)))
        elseif op isa QuantumGate
            layers[i] = (; pos=qubits(op), U=matrix(op), param_idx=Int[], deriv=(k->nothing))
        else
            throw(ArgumentError("state-vector gradient only supports QuantumGate or parameterized gates; got $(typeof(op))"))
        end
    end
    return layers
end

# ─────────────────────────────────────────────────────────────────────────────
# Density-Matrix (noisy) version: represent gates/channels as local superoperators S acting on 2n bits
#   - Gate U  →  S = conj(U) ⊗ U
#   - Map {Kℓ} → S = Σℓ conj(Kℓ) ⊗ Kℓ
#   - Parameter derivatives: use central difference at the superoperator level
# ─────────────────────────────────────────────────────────────────────────────

# Construct gate superoperator
_superop_gate(U::AbstractMatrix) = kron(conj.(U), U)

# Construct channel superoperator（Σℓ conj(Kℓ) ⊗ Kℓ）
function _superop_map(Ks::AbstractVector{<:AbstractMatrix})
    d = size(Ks[1], 1)
    A = zeros(eltype(Ks[1]), d*d, d*d)
    @inbounds for K in Ks
        A .+= kron(conj.(K), K)
    end
    return A
end

# For param-op: use central difference at superoperator level
function _local_derivative_superop(op, k::Int; eps::Float64=1e-6)
    p⁺ = copy(op.params); p⁻ = copy(op.params)
    p⁺[k] += eps; p⁻[k] -= eps
    X⁺ = op.build(p⁺); X⁻ = op.build(p⁻)
    S⁺ = X⁺ isa QuantumGate ? _superop_gate(matrix(X⁺)) : _superop_map(kraus(X⁺))
    S⁻ = X⁻ isa QuantumGate ? _superop_gate(matrix(X⁻)) : _superop_map(kraus(X⁻))
    return (S⁺ .- S⁻) ./ (2eps)
end

# Convert circuit into array of superoperator layers;
# each layer stores (pos2, S, param_idx, deriv)
# Note: the DM version needs total qubit count n to form pos2 = vcat(pos, n .+ pos)
function _materialize_layers_dm(c::QCircuit, n::Int)
    ops = c
    layers = Vector{NamedTuple}(undef, length(ops))
    for (i, op) in enumerate(ops)
        if op isa ParamOp
            X = op.build(op.params)
            if X isa QuantumGate
                pos = qubits(X)
                S   = _superop_gate(matrix(X))
            else
                pos = qubits(X)      
                S   = _superop_map(kraus(X))
            end
            pos2 = vcat(pos, n .+ pos)

            idx = Int[]
            for k in eachindex(op.params)
                op.mask[k] && push!(idx, k)
            end
            layers[i] = (; pos=pos2, S, param_idx=idx,
                         deriv=(k-> _local_derivative_superop(op, k)))
        elseif op isa QuantumGate
            pos2 = vcat(qubits(op), n .+ qubits(op))
            S    = _superop_gate(matrix(op))
            layers[i] = (; pos=pos2, S, param_idx=Int[], deriv=(k->nothing))
        elseif op isa QuantumMap
            pos2 = vcat(qubits(op), n .+ qubits(op))
            S    = _superop_map(kraus(op))
            layers[i] = (; pos=pos2, S, param_idx=Int[], deriv=(k->nothing))
        else
            throw(ArgumentError("density-matrix gradient only supports QuantumGate/QuantumMap/param-op；got $(typeof(op))"))
        end
    end
    return layers
end

# -----------------------------------------------------------------------------
# Two-state backward pass: compute gradient of pure-state expectation
# ∂L/∂θ_j = 2 Re ⟨Ψ_j | (∂Q_j/∂θ_j) | ψ_{j-1}⟩
# where |Ψ_M⟩ = H |ψ_M⟩, and ψ_{j-1} is obtained by inverse evolution (using Q_j^{-1})
# -----------------------------------------------------------------------------

function _gradient_expectation_svec(op::QubitsOperator, circ::QCircuit, ψ0::StateVector)
    layers = _materialize_layers(circ)

    # Forward final state ψ_M
    ψM = StateVector{eltype(ψ0)}(copy(ψ0.data))
    for L in layers
        _apply_gate!(ψM.data, L.pos, L.U)
    end

    # Ψ_M = H ψ_M
    Ψ = _apply_operator_on_state(op, ψM)

    # State to backpropagate (“current ψ”, starting from ψ_M and rolling back)
    Φ = StateVector{eltype(ψ0)}(copy(ψM.data))

    grads = Float64[]   # Order matches active parameter appearance
    for j = length(layers):-1:1
        L = layers[j]

        # Roll back Φ: ψ_{j-1} = Q_j^{-1} ψ_j
        # For noiseless gates, U^{-1} = U'; here we generally use inv(U)
        _apply_gate!(Φ.data, L.pos, inv(L.U))

        # Accumulate contributions for all active parameters in this layer
        for k in L.param_idx
            dU = L.deriv(k)                               # local derivative matrix ∂Q_j/∂θ
            v  = _apply_local_out(Φ, L.pos, dU)           # v = (∂Q_j/∂θ) ψ_{j-1}
            push!(grads, 2 * real(dot(Ψ.data, v.data)))   # 2 Re ⟨Ψ_j | v⟩
        end

        # Roll back Ψ：Ψ_{j-1} = Q_j' Ψ_j
        _apply_gate!(Ψ.data, L.pos, inv(L.U))
    end

    return reverse(grads) 
end

# ∂L/∂θ_j = ⟨⟨Ψ_j | (∂Ŝ_j/∂θ_j) | Φ_{j-1} ⟩⟩
# where |Φ_M⟩⟩ = |ρ_M⟩⟩， |Ψ_M⟩⟩ = Ĥ |I⟩⟩，backprop uses inv(Ŝ_j)
# (note: channels are generally non-unitary)
function _gradient_expectation_dmat(op::QubitsOperator, circ::QCircuit, ρ0::DensityMatrix)
    n = nqubits(ρ0)
    layers = _materialize_layers_dm(circ, n)

    # forward：ρ_M
    vρM = copy(ρ0.data)
    for L in layers
        _apply_gate!(vρM, L.pos, L.S)
    end

    # |Ψ_M⟩⟩ = Ĥ|I⟩⟩
    # Construct |I⟩⟩ (vectorized identity)
    Tρ = eltype(ρ0)
    d  = 1 << n
    vI = zeros(Tρ, d*d)
    @inbounds for k in 1:d
        vI[1 + (k-1)*(d+1)] = one(Tρ)   # diagonal elements
    end

    vΨ = _apply_operator_on_dmvec(op, n, vI)

    # Backprop buffer: Φ starts from ρ_M and rolls back
    vΦ = copy(vρM)

    grads = Float64[]
    for j = length(layers):-1:1
        L = layers[j]

        # Roll back Φ: Φ_{j-1} = Ŝ_j^{-1} Φ_j
        Sinv = inv(L.S)
        _apply_gate!(vΦ, L.pos, Sinv)

        # Accumulate contributions from active parameters
        for k in L.param_idx
            dS = L.deriv(k)                      # superoperator derivative (4^m × 4^m)
            v  = copy(vΦ)
            _apply_gate!(v, L.pos, dS)           # v = (∂Ŝ_j/∂θ) Φ_{j-1}
            push!(grads, real(dot(vΨ, v)))   # <Ψ_j | v>
        end

        # Roll back Ψ：Ψ_{j-1} = adjoint(Ŝ) Ψ_j   
        _apply_gate!(vΨ, L.pos, adjoint(L.S))
    end

    return reverse(grads)
end


# -----------------------------------------------------------------------------
# Integration with Zygote: register custom pullbacks for expectation nodes (θ)
# -----------------------------------------------------------------------------

# Primitive node: given θ (parameter vector), fixed circ/op/ψ,
# returns ⟨ψ| C(θ)† op C(θ) |ψ⟩
function _expect_over_params(θ::AbstractVector,
                             circ::QCircuit,
                             op::QubitsOperator,
                             ψ::StateVector)
    cθ = _with_params(circ, θ)
    return expectation(op, apply(cθ, ψ))
end

# Custom adjoint: pullback(Δ) returns gradients for relevant variables
Zygote.@adjoint function _expect_over_params(θ::AbstractVector,
                                             circ::QCircuit,
                                             op::QubitsOperator,
                                             ψ::StateVector)
    L = _expect_over_params(θ, circ, op, ψ)
    function pullback(Δ)
        # Use manual backward pass for θ; other variables are treated as constants
        gθ = _gradient_expectation_svec(op, _with_params(circ, θ), ψ)
        return (Δ .* gθ, NoTangent(), NoTangent(), NoTangent())
    end
    return L, pullback
end

# Primitive: given θ, fixed circ/op/ρ, returns Tr[ op * (C(θ) ρ C(θ)†) ]
function _expect_over_params_dm(θ::AbstractVector,
                                circ::QCircuit,
                                op::QubitsOperator,
                                ρ::DensityMatrix)
    cθ = _with_params(circ, θ)
    return expectation(op, apply(cθ, ρ))
end

Zygote.@adjoint function _expect_over_params_dm(θ::AbstractVector,
                                                circ::QCircuit,
                                                op::QubitsOperator,
                                                ρ::DensityMatrix)
    L = _expect_over_params_dm(θ, circ, op, ρ)
    function pullback(Δ)
        gθ = _gradient_expectation_dmat(op, _with_params(circ, θ), ρ)
        return (Δ .* gθ, NoTangent(), NoTangent(), NoTangent())
    end
    return L, pullback
end

# -----------------------------------------------------------------------------
# Public API: gradient(loss, circ)
# -----------------------------------------------------------------------------

"""
    gradient(loss::Function, circ::QCircuit) -> Vector{Float64}

Return gradients for all active parameters (mask==true) in the circuit.

This function:
- **Specialized mode**:
  (If `loss` is constructed by `LossExpectationRealSV` or `LossExpectationRealDM`,
   non-Hermitian operators are unsupported — results may be incorrect.)
  When `loss` is built via `LossExpectationRealSV` or `LossExpectationRealDM`,
  the gradient is computed using the hand-written backward algorithm implemented above.
  This mode is more efficient and numerically stable — recommended for use.

The result is a `Vector{Float64}` ordered consistently with `active_parameters(circ)`.
"""
#function gradient(loss::Function, circ::QCircuit)
#    θ0 = active_parameters(circ)
#    isempty(θ0) && return Float64[]
#    f = θ -> loss(_with_params(circ, θ))
#    return Zygote.gradient(f, θ0)[1]
#end

# ─────────────────────────────────────────────────────────────────────────────
# Lightweight callable objects: LossExpectationRealSV / LossExpectationRealDM
# Represent explicit forms like
#     loss(circ) = ⟨ψ| C(θ)† op C(θ) |ψ⟩ (and the version which takes the real part)
# ─────────────────────────────────────────────────────────────────────────────

"""
    LossExpectationRealSV(op::QubitsOperator, ψ::StateVector)

Represents `loss(circ) = real(⟨ψ| C(θ)† op C(θ) |ψ⟩)`.
"""
struct LossExpectationRealSV
    op::QubitsOperator
    ψ::StateVector
end
(L::LossExpectationRealSV)(circ::QCircuit) = real(expectation(L.op, apply(circ, L.ψ)))
_as_theta_fun(L::LossExpectationRealSV, circ::QCircuit) =
    θ -> real(_expect_over_params(θ, circ, L.op, L.ψ))

struct LossExpectationRealDM
    op::QubitsOperator
    ρ::DensityMatrix
end
(L::LossExpectationRealDM)(circ::QCircuit) = real(expectation(L.op, apply(circ, L.ρ)))
_as_theta_fun(L::LossExpectationRealDM, circ::QCircuit) =
    θ -> real(_expect_over_params_dm(θ, circ, L.op, L.ρ))

# ─────────────────────────────────────────────────────────────────────────────
# Specialized overloads of gradient:
# when loss is one of our two types, force the use of the hand-written backward pass.
# Any other loss still uses the existing generic method (Zygote's default AD).
# ─────────────────────────────────────────────────────────────────────────────

function gradient(loss::LossExpectationRealSV, circ::QCircuit)
    θ0 = active_parameters(circ)
    isempty(θ0) && return Float64[]
    fθ = _as_theta_fun(loss, circ)        
    return Zygote.gradient(fθ, θ0)[1]
end

function gradient(loss::LossExpectationRealDM, circ::QCircuit)
    θ0 = active_parameters(circ)
    isempty(θ0) && return Float64[]
    fθ = _as_theta_fun(loss, circ)
    return Zygote.gradient(fθ, θ0)[1]
end


