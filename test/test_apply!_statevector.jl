using Test
using LinearAlgebra
using MyJuliVQC
using MyJuliVQC.Gates

# convenient aliases
matrix = MyJuliVQC.matrix
qubits = MyJuliVQC.qubits

# ========== Utility: build a “full-space operator” using the same qubit-ordering convention as apply! ==========

# Permute the qubit order from the original (1..n) to `new_order` (a length-n permutation).
# Convention: qubit 1 is the global LSB (least significant bit).
function permute_qubits_matrix(n::Int, new_order::AbstractVector{<:Integer})
    @assert sort(collect(new_order)) == collect(1:n)
    d = 1 << n
    P = zeros(ComplexF64, d, d)
    @inbounds for j in 0:(d-1)  # original basis index (0-based)
        # Read the bits in the original order
        # Construct the integer index i_new in the new ordering
        i_new = 0
        @inbounds for i_newpos in 1:n
            q_old = new_order[i_newpos]             # the i_newpos-th position in the new order comes from this old qubit
            bit = (j >> (q_old-1)) & 0x1
            i_new |= (bit << (i_newpos-1))
        end
        # Map column j → row i_new
        P[i_new+1, j+1] = 1
    end
    return P
end

# Embed a local operator U (size 2^m) into the full n-qubit Hilbert space,
# following the qubit ordering specified by `positions`.
function embed_full(n::Int, positions::AbstractVector{<:Integer}, U::AbstractMatrix)
    m = length(positions)
    others = setdiff(collect(1:n), positions)
    # permutation: positions occupy the low bits (local LSB → MSB corresponds directly to the order of positions)
    new_order = vcat(collect(positions), others)  

    P = permute_qubits_matrix(n, new_order)
    d_rest = 1 << (n - m)

    U_big = kron(Matrix{eltype(U)}(I, d_rest, d_rest), Matrix(U))  
    return P' * U_big * P
end

# Multiply all gates in the circuit to obtain the global unitary U.
# Only supports noiseless circuits: QuantumGate (no QuantumMap).
function circuit_full_unitary(n::Int, circ::MyJuliVQC.QCircuit)
    # instantiate parametric gates
    ops = CircuitElement[]
    for x in circ
        if x isa MyJuliVQC.Gates.ParamOp
            push!(ops, x.build(x.params))
        else
            push!(ops, x)
        end
    end
    # accumulate product: gates are applied left-to-right, so we left-multiply the state vector
    d = 1 << n
    Uall = Matrix{ComplexF64}(I, d, d)
    for op in ops
        @assert op isa MyJuliVQC.QuantumGate "StateVector evolution supports only noiseless circuits (QuantumGate only, no QuantumMap)"
        pos = qubits(op)
        Uloc = matrix(op)
        Ufull = embed_full(n, pos, Uloc)
        Uall = Ufull * Uall
    end
    return Uall
end

# ========== Build initial state for n=6: [1,2im,3,4im,...,63,64im], normalized ==========
function make_init_state(n::Int)
    d = 1 << n
    v = Vector{ComplexF64}(undef, d)
    for k in 1:d
        if isodd(k)
            v[k] = k + 0im
        else
            v[k] = 0 + k*im
        end
    end
    v ./= norm(v)
    return MyJuliVQC.StateVector(v)
end

@testset "apply! on StateVector matches full-matrix evolution (n=6)" begin
    n = 6
    ψ0 = make_init_state(n)
    ψ1 = MyJuliVQC.StateVector(copy(ψ0.data))  # make a copy for apply!

    # build circuit
    circ = QCircuit(
        [XGate(1),
        RzGate(5, π/7),
        RyGate(3, π/6),
        CNOTGate(3, 2),
        CRzGate(4, 5, π/8),        
        TOFFOLIGate(1, 2, 3),
        FREDKINGate(3, 6, 5)]
    )

    # 1) apply! (in-place)
    MyJuliVQC.apply!(circ, ψ1)

    # 2) reference full-matrix evolution
    U = circuit_full_unitary(n, circ)
    ϕ = U * ψ0.data

    # compare
    @test isapprox(ψ1.data, ϕ; atol=1e-10, rtol=0)
end

