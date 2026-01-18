using Test
using LinearAlgebra
using MyJuliVQC
using MyJuliVQC.Gates

# --- convenient aliases ---
m(g) = MyJuliVQC.matrix(g)
q(g) = MyJuliVQC.qubits(g)

# --- basic 1-qubit matrices (2×2) ---
I2(T=ComplexF64) = Matrix{T}(I, 2, 2)
XU(T=ComplexF64) = T[0 1; 1 0]
YU(T=ComplexF64) = T[0 -im; im 0]
ZU(T=ComplexF64) = T[1 0; 0 -1]
HU(T=ComplexF64) = (one(T)/sqrt(T(2))) * T[1 1; 1 -1]
SU(T=ComplexF64) = T[1 0; 0 im]
RxU(θ; T=ComplexF64) = T[cos(θ/2) -im*sin(θ/2);
                         -im*sin(θ/2)  cos(θ/2)]

# --- 2-qubit Kronecker embedding ---
kronL(U, T=ComplexF64) = kron(I2(T), U)  # acts on the first qubit of the pair
kronR(U, T=ComplexF64) = kron(U, I2(T))  # acts on the second qubit

# --- absorption utility: given a 2-qubit matrix M and a 1-qubit U, 
#     compute the absorbed matrix according to side = :left / :right
#     and whether it acts on the first or second qubit (:p1 / :p2)
# side = :left  → single-qubit gate is applied after the two-qubit gate (left-multiply)
# side = :right → applied before the two-qubit gate  (right-multiply)
function absorb_matrix(M::AbstractMatrix, U::AbstractMatrix; on::Symbol, side::Symbol)
    T = eltype(M)
    K = on === :p1 ? kronL(U, T) : kronR(U, T)
    return side === :left ? (K * M) : (M * K)
end

@testset "fuse_gates: comprehensive sequence absorption" begin
    θ_fsim = π/5
    ϕ_fsim = π/6
    θ_rx   = π/3

    # Original large circuit (contains adjacent/non-adjacent 1q & 2q gates,
    # including parameterized and special-purpose gates)
    circ = QCircuit([
        CNOTGate(1, 2),                           # 1  two-qubit
        XGate(2),                                 # 2  1q (adjacent; can be absorbed by above CNOT)
        ZGate(2),                                 # 3  1q (non-adjacent to CNOT; cannot be absorbed by it)
        FSIMGate(2, 3, θ_fsim, ϕ_fsim; isparas=[true, true]),  # 4 two-qubit (absorbs Z(2), then RxGate(2))
        RxGate(2, θ_rx; isparas=true),                # 5  1q (adjacent on the right of FSIM)
        CRxGate(1, 2, π/4; isparas=false),            # 6  two-qubit (non-parametric, absorbs right H(2))
        HGate(2),                                 # 7  1q
        SGate(1),                                 # 8  1q (on different qubit than upcoming CZ → not absorbed)
        CZGate(2, 3),                             # 9  two-qubit (absorbs right X(3))
        XGate(3),                                 # 10 1q
        TOFFOLIGate(1, 2, 3),                     # 11 three-qubit (never absorbs)
        YGate(1),                                 # 12 1q (absorbed by final CNOT as a left neighbor)
        CNOTGate(2, 1)                            # 13 two-qubit
    ])

    fused = MyJuliVQC.fuse_gates(circ)

    # Expected structure (from left to right):
    # 1) CNOT(1,2) absorbs its right-neighbor X(2) → (X⊗I) * CNOT
    # 2) Z(2) is absorbed by the FSIM and therefore removed
    # 3) FSIM(2,3) absorbs left Z(2) (right-multiplied), 
    #    then absorbs right RxGate(2) (left-multiplied)
    # 4) CRxGate(1,2,π/4) absorbs the right H(2) (left-multiplied)
    # 5) S(1) is untouched
    # 6) CZ(2,3) absorbs right X(3) (left-multiplied)
    # 7) TOFFOLI unchanged
    # 8) Y(1) is absorbed by the final CNOT and thus removed
    # 9) final CNOT(2,1) absorbs left Y(1) (right-multiplied)
    #
    # So fused should contain: 7 elements (all QuantumGate)
    @test fused isa MyJuliVQC.QCircuit
    @test length(fused) == 7
    @test all(x -> x isa MyJuliVQC.QuantumGate, fused.operations)

    # Check qubit positions
    @test q(fused[1]) == [1, 2]     # CNOT(1,2) absorb X(2)
    @test q(fused[2]) == [2, 3]     # FSIM(2,3) absorb Z(2) and RxGate(2)
    @test q(fused[3]) == [1, 2]     # CRxGate(1,2) absorb H(2)
    @test q(fused[4]) == [1]        # S(1)
    @test q(fused[5]) == [2, 3]     # CZ(2,3) absorb X(3)
    @test q(fused[6]) == [1, 2, 3]  # TOFFOLI(1,2,3)
    @test q(fused[7]) == [2, 1]     # CNOT(2,1) absorb Y(1)

    # Numerically validate matrices
    atol = 1e-10

    # 1) (X ⊗ I) * CNOT(1,2)
    M1_expect = absorb_matrix(m(CNOTGate(1,2)), XU(); on=:p2, side=:left)
    @test isapprox(m(fused[1]), M1_expect; atol=atol, rtol=0)

    # 2) FSIM absorbs Z(2) (right), then Rx(2) (left)
    Mfsim = m(FSIMGate(2,3, θ_fsim, ϕ_fsim; isparas=[false,false]))  
    M2_expect = absorb_matrix(absorb_matrix(Mfsim, ZU(); on=:p1, side=:right),
                              RxU(θ_rx); on=:p1, side=:left)
    @test isapprox(m(fused[2]), M2_expect; atol=atol, rtol=0)

    # 3) (H ⊗ I) * CRxGate(1,2,π/4)
    M3_expect = absorb_matrix(m(CRxGate(1,2, π/4; isparas=false)), HU(); on=:p2, side=:left)
    @test isapprox(m(fused[3]), M3_expect; atol=atol, rtol=0)

    # 4) S(1) unchanged
    @test isapprox(m(fused[4]), m(SGate(1)); atol=atol, rtol=0)

    # 5) (X ⊗ I) * CZ(2,3)
    M5_expect = absorb_matrix(m(CZGate(2,3)), XU(); on=:p2, side=:left)
    @test isapprox(m(fused[5]), M5_expect; atol=atol, rtol=0)

    # 6) TOFFOLI unchanged
    @test isapprox(m(fused[6]), m(TOFFOLIGate(1,2,3)); atol=atol, rtol=0)

    # 7) CNOT(2,1) * (Y(1) ⊗ I)
    M7_expect = absorb_matrix(m(CNOTGate(2,1)), YU(); on=:p2, side=:right)
    @test isapprox(m(fused[7]), M7_expect; atol=atol, rtol=0)
end

@testset "fuse_gates: error handling / invalid inputs" begin
    fg = MyJuliVQC.fuse_gates

    # 1) circuit contains QuantumMap (noise channel) → must throw
    c1 = QCircuit([XGate(1), Depolarizing(1; p=0.1), CNOTGate(1, 2)])
    @test_throws ArgumentError fg(c1)

    # 2) circuit contains invalid elements (skipped QCircuit-level validation) → must throw
#    c2 = QCircuit(Any[XGate(1), 123, ZGate(2)])
#    @test_throws ArgumentError fg(c2)

    # 3) “fake parametric op”: build returns non-QuantumGate → must throw
    #fake_paramop = (; params=[π/7], mask=BitVector([true]), build = p -> 3.14)
    #c3 = QCircuit([XGate(1), fake_paramop])
    #@test_throws ArgumentError fg(c3)

    # 4) Nested circuit contains a QuantumMap (detected after flattening) → must throw
    inner = QCircuit([PhaseDamping(1; γ=0.2)])
    c4 = QCircuit([XGate(1), inner, CNOTGate(1, 2)])
    @test_throws ArgumentError fg(c4)
end

