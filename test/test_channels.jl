using Test
using LinearAlgebra
using MyJuliVQC

# convenient aliases
kraus  = MyJuliVQC.kraus
qubits = MyJuliVQC.qubits
arity  = MyJuliVQC.arity

# small helper
I2(T=ComplexF64) = Matrix{T}(I, 2, 2)

# unified trace-preserving tolerance (consistent with internal implementation,
# QuantumMap uses 1e-5)
const TP_ATOL_F32 = 1e-5
const TP_ATOL_F64 = 1e-5  # also 1e-5 here to match internal checks

@testset "Channels: AmplitudeDamping" begin
    # γ = 0 -> identity channel
    M0 = MyJuliVQC.AmplitudeDamping(1; γ=0.0)
    @test qubits(M0) == [1]
    @test arity(M0) == 1
    @test length(kraus(M0)) == 2
    Ks = kraus(M0)
    @test size(Ks[1]) == (2,2)
    @test size(Ks[2]) == (2,2)
    S = Ks[1]'*Ks[1] + Ks[2]'*Ks[2]
    @test isapprox(S, I2(); atol=TP_ATOL_F64, rtol=0)

    # γ = 1 -> |1⟩ always decays to |0⟩，must still be TP
    M1 = MyJuliVQC.AmplitudeDamping(1; γ=1.0)
    Ks = kraus(M1)
    S = Ks[1]'*Ks[1] + Ks[2]'*Ks[2]
    @test isapprox(S, I2(); atol=TP_ATOL_F64, rtol=0)

    # intermediate value + force Float32
    Mmid = MyJuliVQC.AmplitudeDamping(1; γ=0.3, elty=ComplexF32)
    @test eltype(kraus(Mmid)[1]) == ComplexF32
    @test eltype(kraus(Mmid)[2]) == ComplexF32
    Ks = kraus(Mmid)
    S = Ks[1]'*Ks[1] + Ks[2]'*Ks[2]
    @test isapprox(S, I2(ComplexF32); atol=TP_ATOL_F32, rtol=0)

    # invalid γ
    @test_throws ArgumentError MyJuliVQC.AmplitudeDamping(1; γ=-1e-3)
    @test_throws ArgumentError MyJuliVQC.AmplitudeDamping(1; γ=1.001)

    # invalid position
    @test_throws ArgumentError MyJuliVQC.AmplitudeDamping(0; γ=0.2)
end

@testset "Channels: PhaseDamping" begin
    # γ = 0 -> identity channel
    M0 = MyJuliVQC.PhaseDamping(1; γ=0.0)
    Ks = kraus(M0)
    S = Ks[1]'*Ks[1] + Ks[2]'*Ks[2]
    @test isapprox(S, I2(); atol=TP_ATOL_F64, rtol=0)

    # γ = 1 -> full dephasing, still TP
    M1 = MyJuliVQC.PhaseDamping(1; γ=1.0)
    Ks = kraus(M1)
    S = Ks[1]'*Ks[1] + Ks[2]'*Ks[2]
    @test isapprox(S, I2(); atol=TP_ATOL_F64, rtol=0)

    # intermediate value + Float32
    Mmid = MyJuliVQC.PhaseDamping(1; γ=0.55, elty=ComplexF32)
    @test eltype(kraus(Mmid)[1]) == ComplexF32
    Ks = kraus(Mmid)
    S = Ks[1]'*Ks[1] + Ks[2]'*Ks[2]
    @test isapprox(S, I2(ComplexF32); atol=TP_ATOL_F32, rtol=0)

    # invalid γ
    @test_throws ArgumentError MyJuliVQC.PhaseDamping(1; γ=-1e-6)
    @test_throws ArgumentError MyJuliVQC.PhaseDamping(1; γ=1.001)

    # invalid position
    @test_throws ArgumentError MyJuliVQC.PhaseDamping(0; γ=0.2)
end

@testset "Channels: Depolarizing" begin
    # p = 0 → identity (K0 = I, others = 0), still TP
    M0 = MyJuliVQC.Depolarizing(1; p=0.0)
    @test length(kraus(M0)) == 4
    Ks = kraus(M0)
    S = zero(Ks[1])
    for K in Ks
        S .+= K' * K
    end
    @test isapprox(S, I2(); atol=TP_ATOL_F64, rtol=0)

    # p = 4/3 → K0 = 0, the others each have weight √(1/3)·{X,Y,Z}, still TP
    Mmax = MyJuliVQC.Depolarizing(1; p=4/3)
    Ks = kraus(Mmax)
    S = zero(Ks[1])
    for K in Ks
        S .+= K' * K
    end
    @test isapprox(S, I2(); atol=TP_ATOL_F64, rtol=0)

    # intermediate p + Float32
    Mmid = MyJuliVQC.Depolarizing(1; p=0.4, elty=ComplexF32)
    @test all(eltype(K) == ComplexF32 for K in kraus(Mmid))
    Ks = kraus(Mmid)
    S = zero(Ks[1])
    for K in Ks
        S .+= K' * K
    end
    @test isapprox(S, I2(ComplexF32); atol=TP_ATOL_F32, rtol=0)

    # invalid p (<0 or >4/3)
    @test_throws ArgumentError MyJuliVQC.Depolarizing(1; p=-1e-3)
    @test_throws ArgumentError MyJuliVQC.Depolarizing(1; p=1.5)

    # invalid position
    @test_throws ArgumentError MyJuliVQC.Depolarizing(0; p=0.2)
end
