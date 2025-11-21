using Test
using LinearAlgebra
using MyJuliVQC
using MyJuliVQC: QubitsTerm, QubitsOperator, expectation, StateVector, DensityMatrix

const positions = MyJuliVQC.positions
const oplist    = MyJuliVQC.oplist
const coeff     = MyJuliVQC.coeff

# ----------------------------
# helpers: states & utilities
# ----------------------------

# A fixed deterministic complex vector of dimension 64 (6 qubits), then normalized
function psi1_vec()
    v = Vector{ComplexF64}(undef, 64)
    for k in 1:64
        v[k] = isodd(k) ? (k + 0im) : (0 + k*im)
    end
    v ./= norm(v)
    return v
end

function psi2_vec()
    v = Vector{ComplexF64}(undef, 64)
    for k in 1:64
        a = 2*k - 1
        v[k] = isodd(k) ? (a + 0im) : (0 + a*im)
    end
    v ./= norm(v)
    return v
end

# Kronecker for full L-qubit operator with bit-order convention:
# position 1 (lowest qubit) -> RIGHTMOST factor
# i.e., kron(mats[L], ..., mats[2], mats[1])
function _kron_full(mats::AbstractVector{<:AbstractMatrix})
    @assert !isempty(mats)
    M = Matrix(mats[1])
    @inbounds for k in 2:length(mats)
        M = kron(mats[k], M)
    end
    return M
end

# Materialize a QubitsOperator into a dense 2^L × 2^L matrix (for testing only)
function dense_matrix(op::QubitsOperator, L::Int)
    I2 = ComplexF64[1 0; 0 1]
    H  = zeros(ComplexF64, 1<<L, 1<<L)
    for t in op
        mats = [I2 for _ in 1:L]
        for (p, M) in zip(positions(t), oplist(t))
            mats[p] = ComplexF64.(M)
        end
        H .+= coeff(t) .* _kron_full(mats)
    end
    return H
end

# Heisenberg-1D (open boundary) builder
function heisenberg_1d(L; hz=1.0, J=1.0)
    terms = []
    # on-site Z
    for i in 1:L
        push!(terms, QubitsTerm(i=>"Z"; coeff=hz))
    end
    # nearest-neighbour XX + YY + ZZ
    for i in 1:(L-1)
        push!(terms, QubitsTerm(i=>"X", i+1=>"X"; coeff=J))
        push!(terms, QubitsTerm(i=>"Y", i+1=>"Y"; coeff=J))
        push!(terms, QubitsTerm(i=>"Z", i+1=>"Z"; coeff=J))
    end
    return QubitsOperator(terms)
end

@testset "expectation vs dense-matrix baseline (L=6)" begin
    L  = 6
    ψ1 = psi1_vec()
    ψ2 = psi2_vec()
    ρmix = 0.3*(ψ1*ψ1') + 0.7*(ψ2*ψ2')   # already normalized: tr(ρ)=1

    op = heisenberg_1d(L; hz=1.0, J=2.0)
    H  = dense_matrix(op, L)

    # --- pure states ---
    e1_mine   = expectation(op, StateVector(ψ1))
    e1_dense  = dot(ψ1, H*ψ1)           # scalar = ψ1'*(Hψ1)
    @test isapprox(e1_mine, e1_dense; atol=1e-10, rtol=1e-10)

    e2_mine   = expectation(op, StateVector(ψ2))
    e2_dense  = dot(ψ2, H*ψ2)
    @test isapprox(e2_mine, e2_dense; atol=1e-10, rtol=1e-10)

    # --- mixed state ---
    ρ_vec = vec(ρmix)                          # vectorize in column-major order
    ρ     = DensityMatrix(ρ_vec)               

    eρ_mine  = expectation(op, ρ)
    eρ_dense = tr(H * ρmix)                    # tr(Hρ)
    @test isapprox(eρ_mine, eρ_dense; atol=1e-10, rtol=1e-10)
end
