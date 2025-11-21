using LinearAlgebra
"""
        QuantumGate{T}

Carrier type for quantum gates.
`positions` specifies the qubit indices (1-based) on which the gate acts,
and `data` is the `2^m × 2^m` unitary matrix (column-major), where
`m = length(positions)`.
"""
    struct QuantumGate{T<:Number} <: QuantumPrimitive
        positions::Vector{Int}   # Not sorted; preserves the user-given qubit order
        data::Matrix{T}          # Column-major matrix

        """
        QuantumGate(positions, data; validate=false, atol=1e-10, rtol=0)

        Construct a quantum gate:

        - `positions`: the qubit indices this gate acts on  
        (`Tuple` or `Vector`, 1-based, must be distinct)
        - `data`: a `2^m × 2^m` unitary matrix (column-major),  
        where `m = length(positions)`

        Validation:
        - `positions` must be nonempty, all ≥ 1, and all distinct
        - `data` must be square, its dimension must be a power of 2,  
        and must match `length(positions)`
        - When `validate=true`, additionally check unitarity: `data' * data ≈ I`
        """
        function QuantumGate{T}(positions::AbstractVector{<:Integer},
                            data::AbstractMatrix{T};
                            validate::Bool=false, atol::Real=1e-10, rtol::Real=0
                           ) where {T<:Number}

            # ---- Validate positions ----
            pos = collect(Int, positions)
            isempty(pos) && throw(ArgumentError("positions must be non-empty"))
            any(p -> p < 1, pos) && throw(ArgumentError("qubit indices must be ≥ 1 (1-based)"))
            length(unique(pos)) != length(pos) && throw(ArgumentError("qubit indices must be unique"))
            m = length(pos)

            # ---- Validate data shape ----
            size(data, 1) == size(data, 2) || throw(ArgumentError("gate matrix must be square"))
            d = size(data, 1)
            ispow2(d) || throw(ArgumentError("gate matrix size must be 2^m, got $d"))
            (d == (1 << m)) || throw(ArgumentError("gate matrix size $d does not match positions length m=$m"))

            # Copy into a column-major Matrix and fix element type
            M = Matrix{T}(data)

            # ---- Optional unitarity validation ----
            if validate
                Id = Matrix{T}(LinearAlgebra.I, d, d)
                @assert isapprox(M' * M, Id; atol=atol, rtol=rtol) "gate matrix is not unitary"
            end

            return new{T}(pos, M)
        end
    end

    # Convenience queries
    arity(g::QuantumGate) = length(g.positions)              # Number of qubits m
    matrix(g::QuantumGate) = g.data                          # Extract the matrix
    qubits(g::QuantumGate) = g.positions                     # Extract the qubit index vector

    # Entry for Vector / AbstractVector 
    QuantumGate(positions::AbstractVector{<:Integer},
                data::AbstractMatrix{T}; kwargs...) where {T<:Number} =
        QuantumGate{T}(positions, data; kwargs...)

    # Convenience overload for Tuple positions
    QuantumGate(positions::NTuple{N,<:Integer},
                data::AbstractMatrix{T}; kwargs...) where {N,T<:Number} =
        QuantumGate{T}(collect(positions), data; kwargs...)

    # Convenience overload forcing element type
    QuantumGate(::Type{T}, positions, data; kwargs...) where {T<:Number} =
        QuantumGate{T}(collect(Int, positions), Matrix{T}(data); kwargs...)