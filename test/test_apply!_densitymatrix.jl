using Test
using LinearAlgebra
using MyJuliVQC
using MyJuliVQC.Gates
using MyJuliVQC.Channels

# convenient aliases
matrix = MyJuliVQC.matrix
qubits = MyJuliVQC.qubits
kraus  = MyJuliVQC.kraus
#flattened_ops = MyJuliVQC.flattened_ops

# ---------------- Utility: embed a local operator into the n-qubit full Hilbert space ----------------

# qubit 1 is the global LSB. Given a permutation new_order (length n),
# construct the corresponding 2^n×2^n permutation matrix.
function permute_qubits_matrix(n::Int, new_order::AbstractVector{<:Integer})
    @assert sort(collect(new_order)) == collect(1:n)
    d = 1 << n
    P = zeros(ComplexF64, d, d)
    @inbounds for j in 0:(d-1)   # column index in the original ordering (0-based)
        i_new = 0
        @inbounds for i_newpos in 1:n
            q_old = new_order[i_newpos]
            bit = (j >> (q_old-1)) & 0x1
            i_new |= (bit << (i_newpos-1))
        end
        P[i_new+1, j+1] = 1
    end
    return P
end

# Embed a local operator U (2^m × 2^m) into the full n-qubit space
# based on the ordering defined by `positions`
function embed_full(n::Int, positions::AbstractVector{<:Integer}, U::AbstractMatrix)
    m = length(positions)
    others = setdiff(collect(1:n), positions)
    new_order = vcat(collect(positions), others)    # positions are placed in the global LSB slots
    P = permute_qubits_matrix(n, new_order)
    d_rest = 1 << (n - m)
    U_big = kron(Matrix{eltype(U)}(I, d_rest, d_rest), Matrix(U))
    return P' * U_big * P
end

function full_evolve_density(n::Int, circ::MyJuliVQC.QCircuit, M0::AbstractMatrix)
    M = Matrix{ComplexF64}(M0)
    for op in circ
        # instantiate parametric gates/channels
        x = op
        if  x isa MyJuliVQC.Gates.ParamOp
            x = x.build(x.params)
        end

        if x isa MyJuliVQC.QuantumGate
            pos = qubits(x)
            U   = matrix(x)
            Ufull = embed_full(n, pos, U)
            M = Ufull * M * adjoint(Ufull)
        elseif x isa MyJuliVQC.QuantumMap
            pos = qubits(x)
            Ks  = kraus(x)
            Msum = zeros(eltype(M), size(M,1), size(M,2))
            for K in Ks
                Kfull = embed_full(n, pos, K)
                Msum .+= Kfull * M * adjoint(Kfull)
            end
            M = Msum
        else
            @assert false "unsupported op in circ: $(typeof(x))"
        end
    end
    return M
end

# ---------------- Initial states ----------------

# |ψ1> = normalized [1, 2im, 3, 4im, ..., 15, 16im]  (n=4, length 16)
function psi1_vec()
    v = Vector{ComplexF64}(undef, 16)
    for k in 1:16
        v[k] = isodd(k) ? (k + 0im) : (0 + k*im)
    end
    v ./= norm(v)
    return v
end

# |ψ2> = normalized [1, 3im, 5, 7im, ..., 29, 31im]  (odd sequence, even positions × i)
function psi2_vec()
    v = Vector{ComplexF64}(undef, 16)
    for k in 1:16
        a = 2k - 1                      # 1,3,5,...,31
        v[k] = isodd(k) ? (a + 0im) : (0 + a*im)
    end
    v ./= norm(v)
    return v
end

# ---------------- Circuit ----------------

let
    n = 4

    # Kraus channel：{ 0.5*CNOT, sqrt(3/4)*CZ } acting on (1,2)
    K1 = 0.5 * matrix(CNOTGate(1,2))
    K2 = sqrt(3/4) * matrix(CZGate(1,2))
    chan12 = MyJuliVQC.QuantumMap([1,2], [K1, K2]; check_tp=true)

    circ = QCircuit([
        Depolarizing(4; p=0.3),     # channel
        XGate(3),                   # unitary gate
        CRxGate(4, 3, π/7; isparas=true),      # parametric gate
        CRzGate(3, 4, π/8; isparas=false),     # fixed gate
        TOFFOLIGate(1, 3, 2),       # 3-qubit gate
        chan12                       # custom Kraus channel
    ])

    # -------- case 1: |ψ1><ψ1| --------
    ψ1 = psi1_vec()
    M1 = ψ1 * ψ1'                       # 2^n × 2^n
    ρ1 = MyJuliVQC.DensityMatrix(M1)   
    ρ1_apply = MyJuliVQC.DensityMatrix(copy(ρ1.data))  # make a copy for apply!

    # run apply!
    MyJuliVQC.apply!(circ, ρ1_apply)

    # reference
    M1_ref = full_evolve_density(n, circ, M1)

    @testset "apply!(DensityMatrix) — pure state projector" begin
        @test isapprox(ρ1_apply.data, vec(M1_ref); atol=1e-10, rtol=0)
    end

    # -------- case 2: 0.3|ψ1><ψ1| + 0.7|ψ2><ψ2| --------
    ψ2 = psi2_vec()
    Mmix = 0.3*(ψ1*ψ1') + 0.7*(ψ2*ψ2')
    ρmix = MyJuliVQC.DensityMatrix(Mmix)
    ρmix_apply = MyJuliVQC.DensityMatrix(copy(ρmix.data))

    # apply!
    MyJuliVQC.apply!(circ, ρmix_apply)

    # reference
    Mmix_ref = full_evolve_density(n, circ, Mmix)

    @testset "apply!(DensityMatrix) — classical mixture of two pure states" begin
        @test isapprox(ρmix_apply.data, vec(Mmix_ref); atol=1e-10, rtol=0)
    end
end
