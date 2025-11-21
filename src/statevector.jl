using LinearAlgebra

#### 1) StateVector type ####

    """
        StateVector{T}

    Container type for an n-qubit pure state.  
    The field `data::Vector{T}` must have length `2^n`.
    
    Implementation-wise, the data is stored as a **flattened vector** in
    column-major order (least-significant index varies fastest).  
    Mathematically, it corresponds to a rank-`n` tensor.
    """
    struct StateVector{T<:Number}
        data::Vector{T}          # flattened array of length 2^n
        # Internal constructor: any direct construction must validate length = power of 2
        function StateVector{T}(data::Vector{T}) where {T<:Number}
            L = length(data)
            @assert ispow2(L) "StateVector length must be a power of 2, got $L"
            return new{T}(data)
        end
    end

    # Basic query utilities: state length, number of qubits, element type
    Base.length(ψ::StateVector) = length(ψ.data)
    Base.eltype(ψ::StateVector{T}) where {T} = T

    """
        nqubits(ψ::StateVector) -> Int

    Return the number of qubits `n` in the state vector, requiring that
    `length(ψ) = 2^n`.
    """
    function nqubits(ψ::StateVector)
        L = length(ψ)
        @assert ispow2(L) "StateVector length must be a power of 2, got $L"
        return floor(Int, log2(L))
    end

    
    # Convenient indexing
    Base.getindex(ψ::StateVector, i::Int) = ψ.data[i]
    Base.setindex!(ψ::StateVector, v, i::Int) = (ψ.data[i] = v)

    
    #### 2) Constructors ####

    # (a) Specify n qubits; initialize to |0...0⟩; default scalar type ComplexF64
    """
    StateVector(n::Integer)
    StateVector(T::Type, n::Integer)

    Construct an `n`-qubit pure state initialized to `|0⟩^{⊗n}`.
    - Without a type parameter, the scalar type defaults to `ComplexF64`.
    - Use forms like `StateVector(Float64, n)` to *explicitly* force the element type.
    """
    function StateVector{T}(n::Integer) where {T<:Number}
    @assert n ≥ 0 "number of qubits must be nonnegative"
    data = zeros(T, 2^n)
    data[1] = one(T)                  # first amplitude of |0...0⟩ is 1
    return StateVector{T}(data)
    end
    StateVector(n::Integer) = StateVector{ComplexF64}(n)
    StateVector(::Type{T}, n::Integer) where {T<:Number} = StateVector{T}(n)

    # (b) Construct directly from a given vector (no normalization unless requested)
    """
    StateVector(data::AbstractVector{T}; normalize::Bool=false)
    StateVector(T::Type, data::AbstractVector)

    Construct a pure state directly from the given array.  
    `length(data)` must be a power of 2.

    By default **no normalization** is performed; pass `normalize=true`
    to normalize the vector.
    """
    function StateVector(data::AbstractVector{T}; normalize::Bool=false) where {T<:Number}
    L = length(data)
    @assert ispow2(L) "length of data must be power of 2, got $L"
    v = Vector{T}(data)               # copy into a contiguous Vector
    if normalize
        nrm = LinearAlgebra.norm(v)
        nrm ≈ 0 ? nothing : (v ./= nrm)
    end
    return StateVector{T}(v)
    end
    StateVector(::Type{T}, data::AbstractVector) where {T<:Number} = 
        StateVector(Vector{T}(data))      # convert to the target scalar type
