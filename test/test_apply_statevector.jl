using Test
using LinearAlgebra
using MyJuliVQC
using MyJuliVQC.Gates

# convenient aliases
matrix = MyJuliVQC.matrix
qubits = MyJuliVQC.qubits

# ========== Utility: construct a "full-space operator" consistent with the bit-ordering used in apply! ==========

# Permute qubit order from original (1..n) to new_order (a length-n permutation),
# and construct the corresponding 2^n×2^n permutation matrix P.
# Convention: qubit 1 is the global LSB (least significant bit).
function permute_qubits_matrix(n::Int, new_order::AbstractVector{<:Integer})
    @assert sort(collect(new_order)) == collect(1:n)
    d = 1 << n
    P = zeros(ComplexF64, d, d)
    @inbounds for j in 0:(d-1)  # basis index in the *original* order (0-based)
        # Read each qubit bit from the original index j
        # Construct the index i_new under the permuted bit order
        i_new = 0
        @inbounds for i_newpos in 1:n
            q_old = new_order[i_newpos]             # which original qubit sits in new position i_newpos
            bit = (j >> (q_old-1)) & 0x1
            i_new |= (bit << (i_newpos-1))
        end
        # column j maps to row i_new
        P[i_new+1, j+1] = 1
    end
    return P
end

# Embed a local U (size 2^m) into the n-qubit full space according to `positions`.
function embed_full(n::Int, positions::AbstractVector{<:Integer}, U::AbstractMatrix)
    m = length(positions)
    others = setdiff(collect(1:n), positions)
    # After permutation, `positions` occupy the lower bits in the order given.
    new_order = vcat(collect(positions), others)  

    P = permute_qubits_matrix(n, new_order)
    d_rest = 1 << (n - m)

    U_big = kron(Matrix{eltype(U)}(I, d_rest, d_rest), Matrix(U))  
    return P' * U_big * P
end

# Multiply all gates in a QCircuit to obtain the global unitary U (noise-free only: QuantumGate).
function circuit_full_unitary(n::Int, circ::MyJuliVQC.QCircuit)
    # First materialize parameterized gates into concrete gates
    ops = CircuitElement[]
    for x in circ
        if x isa MyJuliVQC.Gates.ParamOp
            push!(ops, x.build(x.params))
        else
            push!(ops, x)
        end
    end

    # Accumulate the product.
    # (Left-to-right action on the state ⇒ left-multiply column vector in the same order)
    d = 1 << n
    Uall = Matrix{ComplexF64}(I, d, d)
    for op in ops
        @assert op isa MyJuliVQC.QuantumGate "StateVector supports only noiseless circuits (QuantumMap not allowed)"
        pos = qubits(op)
        Uloc = matrix(op)
        Ufull = embed_full(n, pos, Uloc)
        Uall = Ufull * Uall
    end
    return Uall
end

# ========== Construct initial state (n=6): [1,2im,3,4im,...,63,64im] normalized ==========
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

    # Build circuit
    circ = QCircuit(
        [XGate(1),
        RzGate(5, π/7),
        RyGate(3, π/6),
        CNOTGate(3, 2),
        CRzGate(4, 5, π/8),        
        TOFFOLIGate(1, 2, 3),
        FREDKINGate(3, 6, 5)]
    )

    # 1) apply 
    ψ1 = MyJuliVQC.apply(circ, ψ0)
    ψ2 = circ * ψ0  # overloaded *
    ψ3 = MyJuliVQC.apply_ad(circ, ψ0)

    # 2) Full-matrix evolution
    U = circuit_full_unitary(n, circ)
    ϕ = U * ψ0.data

    # Compare results
    @test isapprox(ψ1.data, ϕ; atol=1e-10, rtol=0)
    @test isapprox(ψ2.data, ϕ; atol=1e-10, rtol=0)
    @test isapprox(ψ3.data, ϕ; atol=1e-10, rtol=0)
end

