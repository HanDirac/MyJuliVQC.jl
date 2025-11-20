module Gates

using LinearAlgebra
using ..MyJuliVQC: QuantumGate

# --- Utilities ---------------------------------------------------------------

# 2×2 identity matrix (with specified element type)
_I2(::Type{T}) where {T} = Matrix{T}(I, 2, 2)
# Compute e^{iφ}
_phase(::Type{T}, φ) where {T} = exp(one(T) * im * T(φ))

# Normalize `isparas` into a BitVector having the same shape as the parameters
#_asmask(isparas, n::Int) = isparas isa Bool ? trues(n) .& isparas : BitVector(isparas)

struct ParamOp{N,T,F}
    # N: number of qubits the gate acts on; T: numeric type; F: builder function type
    params::Vector{T}
    mask::BitVector
    build::F   # (p::AbstractVector{T}) -> QuantumGate{N,T}
end
# Make it callable
# count(identity, ...) counts the number of trues in mask
(op::ParamOp)(p::AbstractVector) = (
    @assert length(p) >= count(identity, op.mask) "need $(count(identity, op.mask)) params, got $(length(p))";
    op.build(p)
)

function _paramify(::Type{T}, paras, isparas, builder::Function) where {T}
    p0   = paras isa Tuple ? collect(T, paras) :
           paras isa AbstractVector ? T.(paras) : [T(paras)]
    mask = isparas isa Bool ? fill(isparas, length(p0)) : BitVector(isparas)

    if !any(mask)
        # No mutable/variational parameters: directly construct a fixed gate
        return builder(p0)  # returns QuantumGate{N,T}
    else
        return ParamOp{nothing,T,typeof(builder)}(p0, BitVector(mask), builder)
    end
end

# Generic: lift a single-qubit gate U (2×2) to a controlled gate
# (control first, target second)
function CONTROLGate(ctrl::Int, tgt::Int; U::AbstractMatrix{S}, elty::Type=ComplexF64) where {S}
    T = elty
    mat = T[1 0 0 0;
            0 U[1,1] 0 U[1,2];
            0 0 1 0;
            0 U[2,1] 0 U[2,2]]
    return QuantumGate((ctrl, tgt), mat)
end

# Generic: lift a single-qubit gate U to a double-controlled gate
# (two controls first, target last)
function CONTROLCONTROLGate(c1::Int, c2::Int, tgt::Int; U::AbstractMatrix{S}, elty::Type=ComplexF64) where {S}
    T = elty
    mat = T[
        1 0 0 0 0 0 0 0;
        0 1 0 0 0 0 0 0;
        0 0 1 0 0 0 0 0;
        0 0 0 U[1,1] 0 0 0 U[1,2];
        0 0 0 0 1 0 0 0;
        0 0 0 0 0 1 0 0;
        0 0 0 0 0 0 1 0;
        0 0 0 U[2,1] 0 0 0 U[2,2]
    ]
    return QuantumGate((c1, c2, tgt), mat)
end

# --- Standard single-qubit gates -------------------------------------------

function XGate(n::Int; T::Type=ComplexF64)
    mat = T[0 1; 1 0]
    return QuantumGate((n,), mat)
end

function YGate(n::Int; T::Type=ComplexF64)
    mat = T[0 -im; im 0]
    return QuantumGate((n,), mat)
end

function ZGate(n::Int; T::Type=ComplexF64)
    mat = T[1 0; 0 -1]
    return QuantumGate((n,), mat)
end

function HGate(n::Int; T::Type=ComplexF64)
    mat = (1/sqrt(T(2))) * T[1 1; 1 -1]
    return QuantumGate((n,), mat)
end

function SGate(n::Int; T::Type=ComplexF64)
    mat = T[1 0; 0 im]
    return QuantumGate((n,), mat)
end

function TGate(n::Int; T::Type=ComplexF64)
    mat = T[1 0; 0 _phase(T, π/4)]
    return QuantumGate((n,), mat)
end

function sqrtXGate(n::Int; T::Type=ComplexF64)
    mat = (1/2) * T[1+im 1-im; 1-im 1+im]
    return QuantumGate((n,), mat)
end

function sqrtYGate(n::Int; T::Type=ComplexF64)
    mat = (1/2) * T[1+im -1-im; 1+im 1+im]
    return QuantumGate((n,), mat)
end

# --- Single-qubit rotation gates -------------------------------------------

# RxGate(θ)
function RxGate(n::Int, theta=0.0; T::Type=ComplexF64, isparas=false)
    builder = p -> begin
        θ = T(p[1])
        U = T[cos(θ/2)  -im*sin(θ/2);
             -im*sin(θ/2)  cos(θ/2)]
        QuantumGate((n,), U)   
    end
    return _paramify(T, theta, isparas, builder)
end

# RyGate(θ)
function RyGate(n::Int, theta=0.0; T::Type=ComplexF64, isparas=false)
    builder = p -> begin
        θ = T(p[1])
        U = T[cos(θ/2)  -sin(θ/2);
             sin(θ/2)    cos(θ/2)]
        QuantumGate((n,), U)
    end
    return _paramify(T, theta, isparas, builder)
end

# RzGate(θ)
function RzGate(n::Int, theta=0.0; T::Type=ComplexF64, isparas=false)
    builder = p -> begin
        θ = T(p[1])
        U = T[_phase(T, -θ/2)  0;
              0                _phase(T,  θ/2)]
        QuantumGate((n,), U)
    end
    return _paramify(T, theta, isparas, builder)
end

# --- Controlled rotation gates (optionally parametric) ----------------------

# CRxGate(θ) = |0⟩⟨0|⊗I + |1⟩⟨1|⊗RxGate(θ)
function CRxGate(ctrl::Int, tgt::Int, theta=0.0; T::Type=ComplexF64, isparas=false)
    builder = p -> begin
        θ = T(p[1])
        U = T[
            1 0 0 0;
            0 cos(θ/2) 0 -im*sin(θ/2);
            0 0 1 0;  
            0 -im*sin(θ/2) 0  cos(θ/2)
        ]
        QuantumGate((ctrl, tgt), U)  # Bind positions here
    end
    return _paramify(T, theta, isparas, builder)
end

function CRyGate(ctrl::Int, tgt::Int, theta=0.0; T::Type=ComplexF64, isparas=false)
    builder = p -> begin
        θ = T(p[1])
        U = T[
            1 0 0 0;
            0 cos(θ/2) 0 -sin(θ/2);
            0 0 1 0 ;
            0 sin(θ/2) 0 cos(θ/2)
        ]
        QuantumGate((ctrl, tgt), U)  # Bind positions here
    end
    return _paramify(T, theta, isparas, builder)
end

function CRzGate(ctrl::Int, tgt::Int, theta=0.0; T::Type=ComplexF64, isparas=false)
    builder = p -> begin
        θ = T(p[1])
        U = T[
            1 0 0 0;
            0 _phase(T, -θ/2) 0 0;
            0 0 1 0;
            0 0 0 _phase(T,  θ/2)
        ]
        QuantumGate((ctrl, tgt), U)  # Bind positions here
    end
    return _paramify(T, theta, isparas, builder)
end

# --- Two-qubit gates (non-parametric) --------------------------------------

function SWAPGate(q1::Int, q2::Int; T::Type=ComplexF64)
    mat = T[1 0 0 0;
            0 0 1 0;
            0 1 0 0;
            0 0 0 1]
    return QuantumGate((q1, q2), mat)
end

# iSWAP: |01⟩ ↔ |10⟩ and multiply by i on both swapped states
function iSWAPGate(q1::Int, q2::Int; T::Type=ComplexF64)
    mat = T[1 0 0 0;
            0 0 im 0;
            0 im 0 0;
            0 0 0 1]
    return QuantumGate((q1, q2), mat)
end

# Controlled-NOT (CNOT)
function CNOTGate(q1::Int, q2::Int; T::Type=ComplexF64)
    mat = T[1 0 0 0;
            0 0 0 1;
            0 0 1 0;
            0 1 0 0]
    return QuantumGate((q1, q2), mat)
end

# Controlled-Z
function CZGate(q1::Int, q2::Int; T::Type=ComplexF64)
    mat = T[1 0 0 0;
            0 1 0 0;
            0 0 1 0;
            0 0 0 -1]
    return QuantumGate((q1, q2), mat)
end

# --- Three-qubit gates (non-parametric) -------------------------------------

# Toffoli gate (CCNOT)
function TOFFOLIGate(c1::Int, c2::Int, tgt::Int; T::Type=ComplexF64)
    mat = T[
        1 0 0 0 0 0 0 0;
        0 1 0 0 0 0 0 0;
        0 0 1 0 0 0 0 0;
        0 0 0 0 0 0 0 1;
        0 0 0 0 1 0 0 0;
        0 0 0 0 0 1 0 0;
        0 0 0 0 0 0 1 0;  
        0 0 0 1 0 0 0 0   
    ]
    return QuantumGate((c1, c2, tgt), mat)
end

# Fredkin (controlled-SWAP): a common decomposition: I⊗I⊗|0⟩⟨0| + SWAP⊗|1⟩⟨1|
function FREDKINGate(ctrl::Int, tgt1::Int, tgt2::Int; T::Type=ComplexF64)
    mat = T[
        1 0 0 0 0 0 0 0;
        0 1 0 0 0 0 0 0;
        0 0 1 0 0 0 0 0;
        0 0 0 0 0 1 0 0;
        0 0 0 0 1 0 0 0;
        0 0 0 1 0 0 0 0;
        0 0 0 0 0 0 1 0;  
        0 0 0 0 0 0 0 1   
    ]
    return QuantumGate((ctrl, tgt1, tgt2), mat)
end

# --- FSIM (optionally parametric)----------------------------------------------------
# Since the origin of the five-parameter version in Weiyou Liao's paper is unclear for now,
# we implement the Google Cirq form first.
function FSIMGate(q1::Int, q2::Int, theta=0.0, phi=0.0; T::Type=ComplexF64, isparas=[false, false])
    paras = (theta, phi)
    builder = p -> begin
        θ = T(p[1])
        ϕ = T(p[2])
        a = cos(θ)
        b = -im * sin(θ)
        c = _phase(T, -ϕ)  # e^{-i φ}
        U = T[
            1 0 0 0;
            0 a b 0;
            0 b a 0;
            0 0 0 c
        ]
        QuantumGate((q1, q2), U)  
    end
    return _paramify(T, paras, isparas, builder)
end
end 