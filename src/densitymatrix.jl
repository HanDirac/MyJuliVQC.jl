using LinearAlgebra

#### 1) DensityMatrix type ####

    """
        DensityMatrix{T}

    Carrier type for an n-qubit **mixed state**. Internally stored in **column-major order** (the lowest-index bit varies fastest):  
    the vector `data::Vector{T}` must have length `4^n = 2^(2n)`.

    Convention :
    - Mathematically regarded as a rank-`2n` tensor, each index of size 2;
    - The ket indices are labeled 1..n, the bra indices are labeled n+1..2n;
    - The data can also be viewed as a `2^n × 2^n` density matrix (column-major layout).
    """
    struct DensityMatrix{T<:Number}
        data::Vector{T}      # Vectorized one-dimensional array, length 4^n
    end

    # Query utilities: length, qubit count, element type
    Base.length(rou::DensityMatrix) = length(rou.data)
    Base.eltype(rou::DensityMatrix{T}) where {T} = T

    """
        nqubits(rou::DensityMatrix) -> Int

    Return the number of qubits `n` contained in the mixed state.  
    Requires that `length(rou) = 4^n`.
    """
    function nqubits(rou::DensityMatrix)
        L = length(rou)
        @assert ispow2(L) "DensityMatrix length must be a power of 2, got $L"
        n2 = floor(Int, log2(L))              # Expect n2 = 2n
        @assert iseven(n2) "vector length must be 2^(2n), got 2^$n2"
        return n2 ÷ 2
    end

    """
        matview(rou::DensityMatrix) -> AbstractMatrix

    以 `2^n * 2^n` 的矩阵视图查看（共享底层存储，不复制）。
    """
    function matview(rou::DensityMatrix)
        d = 1 << nqubits(rou)                   # 2^n
        return reshape(rou.data, d, d)
    end

    #### 2) Constructors ####

    # (a) Specify n qubits, initialized to |0…0⟩⟨0…0|; default scalar type ComplexF64
    """
    DensityMatrix(n::Integer)
    DensityMatrix(T::Type, n::Integer)

    Construct an `n`-qubit mixed state initialized as `|0⟩⟨0|^{⊗n}`.
    - Without a type argument, the scalar type defaults to `ComplexF64`;
    - You can also explicitly enforce a scalar type, e.g. `DensityMatrix(Float32, n)`.
    """
    function DensityMatrix{T}(n::Integer) where {T<:Number}
        @assert n ≥ 0 "number of qubits must be nonnegative"
        d = 1 << n                             # 2^n
        M = zeros(T, d, d)
        M[1, 1] = one(T)
        return DensityMatrix{T}(vec(M))        # Column-major vectorization
    end
    DensityMatrix(n::Integer) = DensityMatrix{ComplexF64}(n)
    DensityMatrix(::Type{T}, n::Integer) where {T<:Number} = DensityMatrix{T}(n)

    # (b) Construct directly from given data (vector / matrix) with optional validation and normalization
    """
    DensityMatrix(data::AbstractVector{T};
                  validate::Bool=false, normalize_trace::Bool=false)
    DensityMatrix(data::AbstractMatrix{T};
                  validate::Bool=false, normalize_trace::Bool=false)
    DensityMatrix(T::Type, data::AbstractVector)
    DensityMatrix(T::Type, data::AbstractMatrix)

    - If a 1D vector is provided, its length must be `4^n`;
    - If a 2D matrix is provided, its size must be `d×d` with `d = 2^n`.

    - When passing a vector, the matrix is assumed to be input **column by column**;
    - When passing a matrix, input the whole matrix directly.

    Keyword arguments:
    - `normalize_trace=true`: rescale trace to 1;
    - `validate=true`: perform basic density-matrix checks (Hermitian, positive semidefinite, `tr ≈ 1`).
    """
    function DensityMatrix(data::AbstractVector{T};
                           validate::Bool=false, normalize_trace::Bool=false) where {T<:Number}
        L = length(data)
        @assert ispow2(L) "length of vectorized density must be power of 2, got $L"
        n2 = floor(Int, log2(L))
        @assert iseven(n2) "vector length must be 2^(2n), got 2^$n2"
        v = Vector{T}(data)
        ## (Noted on 2025-08-19) Currently, when the input vector length is not 2^(2n),
        ## no exception is thrown — this is a known bug for now.
        if normalize_trace || validate
            d = 1 << (n2 ÷ 2)
            M = reshape(v, d, d)
            if normalize_trace
                tr = real(LinearAlgebra.tr(M))
                tr ≈ 0 ? nothing : (M ./= tr)
            end
            if validate
                @assert isapprox(M, M'; atol=1e-10, rtol=0) "density not Hermitian"
                vals = LinearAlgebra.eigvals(LinearAlgebra.Hermitian(M))
                @assert minimum(real(vals)) ≥ -1e-12
                @assert isapprox(real(LinearAlgebra.tr(M)), 1.0; atol=1e-10)
            end
            v = vec(M)                         # Possibly rescaled
        end
        return DensityMatrix{T}(v)
    end

    function DensityMatrix(data::AbstractMatrix{T};
                           validate::Bool=false, normalize_trace::Bool=false) where {T<:Number}
        size(data, 1) == size(data, 2) || throw(ArgumentError("matrix must be square"))
        d = size(data, 1)
        if !ispow2(d)
            throw(ArgumentError("matrix size must be 2^n, got $d"))
        end
        M = Matrix{T}(data)

        if normalize_trace
            tr = real(LinearAlgebra.tr(M))
            tr ≈ 0 ? nothing : (M ./= tr)
        end
        if validate
            @assert isapprox(M, M'; atol=1e-10, rtol=0) "density not Hermitian"
            vals = LinearAlgebra.eigvals(LinearAlgebra.Hermitian(M))
            @assert minimum(real(vals)) ≥ -1e-12
            @assert isapprox(real(LinearAlgebra.tr(M)), 1.0; atol=1e-10)
        end
        return DensityMatrix{T}(vec(M))
    end

    DensityMatrix(::Type{T}, data::AbstractVector) where {T<:Number} =
        DensityMatrix(Vector{T}(data); validate=false, normalize_trace=false)
    DensityMatrix(::Type{T}, data::AbstractMatrix) where {T<:Number} =
        DensityMatrix(Matrix{T}(data); validate=false, normalize_trace=false)