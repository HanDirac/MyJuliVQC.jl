using Test
using LinearAlgebra
using MyJuliVQC
using MyJuliVQC.Gates
using MyJuliVQC.Channels

# Convenient aliases
matrix = MyJuliVQC.matrix
qubits = MyJuliVQC.qubits
kraus  = MyJuliVQC.kraus
#flattened_ops = MyJuliVQC.flattened_ops

# ---------------- Utility: embed local operators into the full n-qubit Hilbert space ----------------

# Given a permutation new_order (a permutation of 1..n), construct the 2^n×2^n permutation matrix.
function permute_qubits_matrix(n::Int, new_order::AbstractVector{<:Integer})
    @assert sort(collect(new_order)) == collect(1:n)
    d = 1 << n
    P = zeros(ComplexF64, d, d)
    @inbounds for j in 0:(d-1)   # column index in the *original* order (0-based)
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

# Embed a local operator U (2^m×2^m) at `positions` into the full n-qubit Hilbert space.
function embed_full(n::Int, positions::AbstractVector{<:Integer}, U::AbstractMatrix)
    m = length(positions)
    others = setdiff(collect(1:n), positions)
    new_order = vcat(collect(positions), others)    # put `positions` in the global low bits
    P = permute_qubits_matrix(n, new_order)
    d_rest = 1 << (n - m)
    U_big = kron(Matrix{eltype(U)}(I, d_rest, d_rest), Matrix(U))
    return P' * U_big * P
end

function full_evolve_density(n::Int, circ::MyJuliVQC.QCircuit, M0::AbstractMatrix)
    M = Matrix{ComplexF64}(M0)
    for op in circ
        # Materialize parameterized gates / channels
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

# |ψ2> = normalized [1, 3im, 5, 7im, ..., 29, 31im]  
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

    # Kraus channel on qubits (1,2): { 0.5*CNOT, sqrt(3/4)*CZ }
    K1 = 0.5 * matrix(CNOTGate(1,2))
    K2 = sqrt(3/4) * matrix(CZGate(1,2))
    chan12 = MyJuliVQC.QuantumMap([1,2], [K1, K2]; check_tp=true)

    circ = QCircuit([
        Depolarizing(4; p=0.3),     # channel
        XGate(3),                   # gate
        CRxGate(4, 3, π/7; isparas=true),      # parameterized gate
        CRzGate(3, 4, π/8; isparas=false),     # fixed gate (still treated as a gate)
        TOFFOLIGate(1, 3, 2),       # 3-qubit gate
        chan12                       # custom Kraus channel
    ])

    # -------- case 1: |ψ1><ψ1| --------
    ψ1 = psi1_vec()
    M1 = ψ1 * ψ1'                       # 2^n × 2^n
    ρ1 = MyJuliVQC.DensityMatrix(M1)    # stored in vectorized form

    # apply
    ρ1_ = apply(circ, ρ1)
    ρ1__ = circ * ρ1  
    ρ1_ad = MyJuliVQC.apply_ad(circ, ρ1)


    # Full-matrix reference evolution
    M1_ref = full_evolve_density(n, circ, M1)

    @testset "apply(DensityMatrix) — pure state projector" begin
        @test isapprox(ρ1_.data, vec(M1_ref); atol=1e-10, rtol=0)
        @test isapprox(ρ1__.data, vec(M1_ref); atol=1e-10, rtol=0)
        @test isapprox(ρ1_ad.data, vec(M1_ref); atol=1e-10, rtol=0)
    end

    # -------- case 2: 0.3|ψ1><ψ1| + 0.7|ψ2><ψ2| --------
    ψ2 = psi2_vec()
    Mmix = 0.3*(ψ1*ψ1') + 0.7*(ψ2*ψ2')
    ρmix = MyJuliVQC.DensityMatrix(Mmix)

    # apply
    ρmix_ = apply(circ, ρmix)
    ρmix__ = circ * ρmix  
    ρmix_ad = MyJuliVQC.apply_ad(circ, ρmix)

    # reference
    Mmix_ref = full_evolve_density(n, circ, Mmix)

    @testset "apply(DensityMatrix) — classical mixture of two pure states" begin
        @test isapprox(ρmix_.data, vec(Mmix_ref); atol=1e-10, rtol=0)
        @test isapprox(ρmix__.data, vec(Mmix_ref); atol=1e-10, rtol=0)
        @test isapprox(ρmix_ad.data, vec(Mmix_ref); atol=1e-10, rtol=0)
    end
end
