using LinearAlgebra

"""
        QuantumMap{T}

Carrier type for quantum channels.
`positions` specifies the qubit indices being acted on (1-based, all distinct).
`kraus` is a list of Kraus matrices, each of size `2^m × 2^m`,
where `m = length(positions)`.
"""
    struct QuantumMap{T<:Number} <: QuantumPrimitive
        positions::Vector{Int}
        kraus::Vector{Matrix{T}}
        """
        QuantumMap(positions, kraus; elty=ComplexF64, check_tp=true)

    Construct a quantum channel.

    Validation:
    - Each Kraus operator must be `2^m × 2^m`;
    - When `check_tp=true`, verify trace preservation: `∑ K'K ≈ I`.
    """
        function QuantumMap(positions, kraus;
                            elty::Type=ComplexF64,
                            check_tp::Bool=true)

            pos = collect(Int, positions)
            isempty(pos)             && throw(ArgumentError("positions must be non-empty"))
            any(p->p<1, pos)         && throw(ArgumentError("qubit indices must be ≥ 1 (1-based)"))
            length(unique(pos))!=length(pos) && throw(ArgumentError("qubit indices must be unique"))
            m = length(pos); d = 1 << m

            Ks = [Matrix{elty}(K) for K in kraus]
            isempty(Ks) && throw(ArgumentError("kraus list must be non-empty"))

            for (i, K) in enumerate(Ks)
                size(K,1)==size(K,2) || throw(ArgumentError("Kraus #$i is not square"))
                size(K,1)==d         || throw(ArgumentError("Kraus #$i must be $d×$d for m=$m"))
            end

            if check_tp
                S = zeros(elty, d, d)
                for K in Ks
                    S .+= K' * K
                end
                @assert isapprox(S, Matrix{elty}(I, d, d); atol=1e-5, rtol=0) "∑ K'K ≠ I within tolerance"
            end

            return new{elty}(pos, Ks)
        end
    end

    # Convenience queries
    arity(M::QuantumMap)  = length(M.positions)
    kraus(M::QuantumMap)  = M.kraus
    qubits(M::QuantumMap) = M.positions

    # ================= QuantumMap convenience constructors =================
    # 1) positions as NTuple → convert to Vector first
    QuantumMap(positions::NTuple{N,<:Integer},
            kraus; kwargs...) where {N} =
        QuantumMap(collect(positions), kraus; kwargs...)

    # 2) Single-qubit convenience: pass an Int (e.g. pos=1)
    QuantumMap(pos::Integer, kraus; kwargs...) =
        QuantumMap([Int(pos)], kraus; kwargs...)

    # 3) Forcing scalar type: ensure Kraus matrices are of element type T
    QuantumMap(::Type{T},
            positions,
            kraus; kwargs...) where {T<:Number} =
        QuantumMap(positions, kraus; elty=T, kwargs...)