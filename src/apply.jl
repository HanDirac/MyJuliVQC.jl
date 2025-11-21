using LinearAlgebra
using Base.Threads

module _ThreadingUtil
    using LinearAlgebra
    const BLAS_mod = LinearAlgebra.BLAS
    # Utility functions for saving/restoring BLAS thread counts
    set_blas_threads!(n::Int) = try; BLAS_mod.set_num_threads(n); catch; end
    get_blas_threads() = try; BLAS_mod.get_num_threads(); catch; 1; end
end

# Top-level short aliases (declared in a legal position)
const get_blas_threads  = _ThreadingUtil.get_blas_threads
const set_blas_threads! = _ThreadingUtil.set_blas_threads!

#### a) apply! for StateVector (noiseless Schrödinger evolution) ####
    """
        apply!(circ::QCircuit, ψ::StateVector)

    Apply the quantum circuit `circ` to the pure state `ψ` (**in place**).
    Only supports **noiseless** circuits: if the circuit contains a `QuantumMap`, an error is thrown.
    For each parametric operation `(; params, mask, build)` in `circ`,
    the gate is instantiated at runtime by calling `build(params)` before application.
    """
    function apply!(circ::QCircuit, ψ::StateVector)
        # Flatten the circuit and apply gates sequentially
#        ops = flattened_ops(circ)
        ops = circ

        # Type alias
        Tψ = eltype(ψ)

        # Execute each operation
        for (k, op) in enumerate(ops)
            if op isa QuantumMap
                throw(ArgumentError("apply!(StateVector) only supports noiseless circuits (found QuantumMap at position $k)"))
            end
            gate =
                if op isa ParamOp
                    op.build(op.params)
                elseif op isa QuantumGate
                    op
                else
                    throw(ArgumentError("unsupported op at position $k: $(typeof(op)); expected QuantumGate or ParamOp"))
                end

            _apply_gate!(ψ.data, qubits(gate), matrix(gate))
        end

        return ψ
    end

    # ---------- Internal utility: apply an arbitrary m-qubit gate to a state vector ----------

    # Apply an m-qubit unitary matrix U to the qubits specified by `positions`.
    # data: state vector of length 2^n (modified in place)
    function _apply_gate!(data::Vector{Tψ}, positions::Vector{Int}, U::AbstractMatrix{S}) where {Tψ,S}
        m = length(positions)
        @assert m ≥ 1 "gate arity must be ≥ 1"
        N = length(data)
        @assert ispow2(N) "state vector length must be a power of 2, got $N"
        n = floor(Int, log2(N))
        @assert all(1 .≤ positions .≤ n) "gate positions out of range 1..$n"

        # ---- Precomputation ----
        # Strides for target bits (determine jumps within each 2^m local block)
        strides = map(p -> 1 << (p - 1), positions)
        dloc    = 1 << m                      # Local dimension 2^m
        drest   = 1 << (n - m)                # Number of enumerations for remaining bits

        # Convert U to the element type
        Uloc = Matrix{Tψ}(U)

        # Precompute global offset off(t) for each local index t ∈ 0..2^m-1
        offs = Vector{UInt}(undef, dloc)
        @inbounds for t in 0:(dloc-1)
            off = UInt(0); tt = t
            @inbounds for k = 1:m
                if (tt & 0x1) == 0x1
                    off += UInt(strides[k])
                end
                tt >>= 1
            end
            offs[t+1] = off
        end

        # Remaining bits and their strides (determine base of each “block”)
        others = setdiff(collect(1:n), positions)       # ascending order
        strides_other = map(p -> 1 << (p - 1), others)

        # Number of columns per panel
        panel = min(8, drest)                           

        # Preallocate panel buffers
        VB   = Matrix{Tψ}(undef, dloc, panel)           # Collected local vectors (columns)
        OUTB = similar(VB)                               # Result after multiplying by U

        # Helper: compute base address for block r (target bits = 0)
        base_from_r(r::Int) = begin
            b = UInt(0); rr = r
            @inbounds for k = 1:length(strides_other)
                if (rr & 0x1) == 0x1
                    b += UInt(strides_other[k])
                end
                rr >>= 1
            end
            b
        end

        # ---- Outer loop over blocks (panels) of remaining bits ----
        # Inner loop performs dloc×B small matrix multiplications
        r0 = 0
        while r0 < drest
            B = min(panel, drest - r0)                   # Actual number of columns this round
            use_outer_threads_local = MyJuliVQC.use_outer_threads() &&
                          (dloc <= MyJuliVQC.dloc_thread_threshold()) &&
                          (B >= 2)
            begin
                # gather: collect B local vectors into the first B columns of VB
                if use_outer_threads_local
                    Threads.@threads for j = 0:(B-1)
                        base = base_from_r(r0 + j)
                        @inbounds @simd for t = 0:(dloc-1)
                            idx = Int(base + offs[t+1]) + 1
                            VB[t+1, j+1] = data[idx]
                        end
                    end
                else
                    for j = 0:(B-1)
                        base = base_from_r(r0 + j)
                        @inbounds @simd for t = 0:(dloc-1)
                            idx = Int(base + offs[t+1]) + 1
                            VB[t+1, j+1] = data[idx]
                        end
                    end
                end

                # OUTB[:,1:B] = Uloc * VB[:,1:B]

                # Simple heuristic: small U + multiple columns → parallel over j; otherwise rely on BLAS
                old_blas = get_blas_threads()
                if use_outer_threads_local
                    set_blas_threads!(1)  # avoid competition
                end

                @views mul!(OUTB[:, 1:B], Uloc, VB[:, 1:B])

                if use_outer_threads_local
                    set_blas_threads!(old_blas)
                end

                # scatter: write back
                if use_outer_threads_local
                    Threads.@threads for j = 0:(B-1)
                        base = base_from_r(r0 + j)
                        @inbounds @simd for t = 0:(dloc-1)
                            idx = Int(base + offs[t+1]) + 1
                            data[idx] = OUTB[t+1, j+1]
                        end
                    end
                else
                    for j = 0:(B-1)
                        base = base_from_r(r0 + j)
                        @inbounds @simd for t = 0:(dloc-1)
                            idx = Int(base + offs[t+1]) + 1
                            data[idx] = OUTB[t+1, j+1]
                        end
                    end
                end
            r0 += B
            end
        end
        return nothing
    end

#### b) apply! for DensityMatrix (noisy / noiseless) ####
"""
    apply!(circ::QCircuit, ρ::DensityMatrix)

Applies the quantum circuit `circ` to the mixed state `ρ` (**in place**). 
Supports circuits containing both noiseless gates (`QuantumGate`) and noisy channels (`QuantumMap`).  
`ρ` is treated as a pure state of `2n` qubits,  
where the ket indices are `1..n` and bra indices are `n+1..2n` (ket on lower bits, bra on higher bits).
"""
function apply!(circ::QCircuit, ρ::DensityMatrix)
#    ops = flattened_ops(circ)
    ops = circ
    n = nqubits(ρ)
    Tρ = eltype(ρ)

    for (k, op) in enumerate(ops)
        # Parametric operation: instantiate first
        x = op isa ParamOp ? op.build(op.params) : op

        if x isa QuantumGate
            _apply_gate_dm!(ρ.data, n, qubits(x), matrix(x))
        elseif x isa QuantumMap
            _apply_map_dm!(ρ.data, n, qubits(x), kraus(x))
        else
            throw(ArgumentError(
                "apply!(DensityMatrix) expects QuantumGate/QuantumMap/ParamOp; got $(typeof(x)) at position $k"))
        end
    end
    return ρ
end

# -- Internal: apply an m-qubit quantum gate U to a density matrix --
# After vectorization, this becomes a 2m-qubit gate conj(U) ⊗ U.
# data: vec(ρ) (length 4^n)
# n: total number of qubits
# pos: physical positions of the gate (1-based)
# U: local 2^m×2^m matrix
function _apply_gate_dm!(data::Vector{Tρ}, n::Int, pos::Vector{Int}, U::AbstractMatrix{S}) where {Tρ,S}
    m = length(pos)
    # In the vectorized 2n-qubit system, ket occupies [pos] (lower bits), bra at [n .+ pos] (higher bits).
    pos2 = vcat(pos, n .+ pos)
    # Local operator: vec(U ρ U†) = (conj(U) ⊗ U) * vec(ρ)
    Uket = Matrix{Tρ}(U)
    Ubra = conj.(Uket)
    U2   = kron(Ubra, Uket)   # Note order: bra (slow index) first, then ket (fast index)
    _apply_gate!(data, pos2, U2)
    return nothing
end

# -- Internal: apply Kraus channel {Kₗ} to a density matrix --
# After vectorization: ∑ₗ conj(Kₗ) ⊗ Kₗ
function _apply_map_dm!(data::Vector{Tρ}, n::Int, pos::Vector{Int}, Ks::Vector{<:AbstractMatrix}) where {Tρ}
    m   = length(pos)
    d   = 1 << m              # 2^m
    pos2 = vcat(pos, n .+ pos)

    # Aggregate into a 4^m×4^m local matrix: A = ∑ₗ conj(Kₗ) ⊗ Kₗ
    A = zeros(Tρ, d*d, d*d)
    @inbounds for K in Ks
        Kt = Matrix{Tρ}(K)
        A .+= kron(conj.(Kt), Kt)
    end
    _apply_gate!(data, pos2, A)
    return nothing
end

#### c) apply (out-of-place) & circ * state ####

"""
    apply(circ::QCircuit, ψ::StateVector) -> StateVector

Out-of-place version: returns a new pure state after applying `circ`  
without modifying the input `ψ`.  
Equivalent to `circ * ψ`.
"""
function apply(circ::QCircuit, ψ::StateVector)
    ψcopy = StateVector{eltype(ψ)}(copy(ψ.data))  # Deep copy data
    apply!(circ, ψcopy)
    return ψcopy
end

"""
    apply(circ::QCircuit, ρ::DensityMatrix) -> DensityMatrix

Out-of-place version: returns a new density matrix after applying `circ`  
without modifying the input `ρ`.  
Equivalent to `circ * ρ`.
"""
function apply(circ::QCircuit, ρ::DensityMatrix)
    ρcopy = DensityMatrix{eltype(ρ)}(copy(ρ.data))  # Deep copy data (preserving vectorized storage)
    apply!(circ, ρcopy)
    return ρcopy
end

# Operator overloading: circ * state ≡ apply(circ, state)
Base.:*(circ::QCircuit, ψ::StateVector)   = apply(circ, ψ)
Base.:*(circ::QCircuit, ρ::DensityMatrix) = apply(circ, ρ)

# === d) AD-friendly apply (no in-place writes visible to AD) ===
using Zygote: Buffer, copy, ignore

"""
    apply_ad(circ::QCircuit, ψ::StateVector) -> StateVector

Same semantics as `apply(circ, ψ)`, but Zygote-friendly:  
internally performs no in-place modifications visible to automatic differentiation.  
Suitable for general AD use cases; slightly slower than `apply!/apply`.
"""

function _apply_gate_ad(data::Vector{T}, positions::Vector{Int}, U::AbstractMatrix{S}) where {T,S}
    m = length(positions)
    N = length(data)
    n = floor(Int, log2(N))
    strides = map(p -> 1 << (p - 1), positions)
    dloc = 1 << m
    drest = 1 << (n - m)

    Uloc = Matrix{T}(U)

    offs = ignore() do
        # Purely functional construction, no setindex!; the entire result is ignored by AD (treated as constant)
        [UInt(sum(((t >> (k-1)) & 0x1) == 1 ? strides[k] : 0 for k in 1:m)) for t in 0:(dloc-1)]
    end

    others = ignore() do
    setdiff(collect(1:n), positions)
    end
    strides_other = map(p -> 1 << (p - 1), others)

    base_from_r_ad(rr, strides_other) = begin
        b = UInt(0); r = rr
        @inbounds for k = 1:length(strides_other)
            if (r & 0x1) == 0x1
                b += UInt(strides_other[k])
            end
            r >>= 1
        end
        b
    end

    panel = min(8, drest)
    next = Zygote.Buffer(similar(data))   # Keep this outside the while loop; copy only once at the end of the layer

    r0 = 0
    while r0 < drest
        B = min(panel, drest - r0)

        # For each panel, create a new Buffer and discard it after use (it will not be written in the next round)
        VBbuf = Zygote.Buffer(Matrix{T}(undef, dloc, B))

        # gather -> VBbuf  (parallelizable: columns of VBbuf accessed by different j do not conflict)
        # Decide whether to enable outer-level parallelism (use a local variable to avoid treating a function name as a Bool)
        use_outer_threads_local =
            MyJuliVQC.use_outer_threads() &&
            (dloc <= MyJuliVQC.dloc_thread_threshold()) &&
            (B >= 2)
        if use_outer_threads_local 
            Threads.@threads for j = 0:(B-1)
                base = base_from_r_ad(r0 + j, strides_other)
                @inbounds @simd for t = 0:(dloc-1)
                    idx = Int(base + offs[t+1]) + 1
                    VBbuf[t+1, j+1] = data[idx]
                end
            end
        else
            for j = 0:(B-1)
                base = base_from_r_ad(r0 + j, strides_other)
                @inbounds @simd for t = 0:(dloc-1)
                    idx = Int(base + offs[t+1]) + 1
                    VBbuf[t+1, j+1] = data[idx]
                end
            end
        end
        # Materialize once before multiplication; VBbuf is frozen immediately, but we no longer write to it this round
        VB = copy(VBbuf)                  # Matrix{T} (dloc, B)

        # Simple heuristic: for small U + many columns, parallelize over j; otherwise delegate to BLAS
        old_blas = get_blas_threads()
        if use_outer_threads_local
            set_blas_threads!(1)  # avoid contention
        end

        OUTsub = Uloc * VB # (dloc, B)

        if use_outer_threads_local
            set_blas_threads!(old_blas)
        end

        # scatter -> Buffer next  (parallelizable: index sets written by different j into next do not overlap)
        if use_outer_threads_local 
            Threads.@threads for j = 0:(B-1)
                base = base_from_r_ad(r0 + j, strides_other)
                @inbounds @simd for t = 0:(dloc-1)
                    idx = Int(base + offs[t+1]) + 1
                    next[idx] = OUTsub[t+1, j+1]
                end
            end
        else
            for j = 0:(B-1)
                base = base_from_r_ad(r0 + j, strides_other)
                @inbounds @simd for t = 0:(dloc-1)
                    idx = Int(base + offs[t+1]) + 1
                    next[idx] = OUTsub[t+1, j+1]
                end
            end
        end
        r0 += B
    end

    return copy(next)  # Only here do we materialize next into an actual Vector
end

function apply_ad(circ::QCircuit, ψ::StateVector)
    # Flattened sequence of operations; keep the same execution order and instantiation logic as apply!
#    ops = flattened_ops(circ)
    ops = circ
    Tψ  = eltype(ψ)

    cur = copy(ψ.data) 

    # Layer by layer: each layer's gather/multiply/scatter writes to a Buffer
    for op in ops
        op isa QuantumMap && throw(ArgumentError(
            "apply_ad(StateVector) only supports noiseless circuits (found QuantumMap at position $k)"))

        gate = op isa ParamOp ? op.build(op.params) :
               (op isa QuantumGate ? op :
                throw(ArgumentError("unsupported operation $(typeof(op)) (expected QuantumGate or ParamOp)")))

        # — Same block algorithm as _apply_gate!, but the target container is a Buffer —
        positions = qubits(gate)
        U   = matrix(gate)
        cur = _apply_gate_ad(cur, positions, U)   # update cur
    end

    return StateVector{Tψ}(copy(cur))  # Materialize into a real Vector and return
end

"""
    apply_ad(circ::QCircuit, ρ::DensityMatrix) -> DensityMatrix

Same semantics as `apply(circ, ρ)`, but Zygote-friendly:  
no in-place modifications are made to arrays visible to differentiation.  
Suitable for “general AD” scenarios (with or without noise).
"""
function apply_ad(circ::QCircuit, ρ::DensityMatrix)
#    ops = flattened_ops(circ)
    ops = circ
    n   = nqubits(ρ)
    Tρ  = eltype(ρ)

    cur = copy(ρ.data)   

    for op in ops
        x = op isa ParamOp ? op.build(op.params) : op

        if x isa QuantumGate
            # Superoperator: conj(U) ⊗ U, acting on the 2n bits at [pos ; n .+ pos]
            pos  = qubits(x)
            Uket = Matrix{Tρ}(matrix(x))
            Ubra = conj.(Uket)
            S    = kron(Ubra, Uket)
            pos2 = vcat(pos, n .+ pos)
            cur = _apply_gate_ad(cur, pos2, S) 
        elseif x isa QuantumMap
            # Aggregate Kraus operators into a superoperator: A = Σₗ conj(Kₗ) ⊗ Kₗ
            pos = qubits(x)
            Ks  = kraus(x)
            m   = length(pos)
            d   = 1 << m

            A = zeros(Tρ, d*d, d*d)
            @inbounds for K in Ks
                Kt = Matrix{Tρ}(K)
                A .+= kron(conj.(Kt), Kt)
            end
            pos2 = vcat(pos, n .+ pos)

            let S = A
                cur = _apply_gate_ad(cur, pos2, S) 
            end
        else
            throw(ArgumentError(
                "apply_ad(DensityMatrix) expected QuantumGate/QuantumMap/ParamOp; got $(typeof(x))"))
        end
    end

    return DensityMatrix{Tρ}(copy(cur))
end

# --- Tell Zygote to ignore in-place kernels (no AD through them) ---
using Zygote

Zygote.@nograd apply!
Zygote.@nograd _apply_gate!
Zygote.@nograd _apply_gate_dm!
Zygote.@nograd _apply_map_dm!

