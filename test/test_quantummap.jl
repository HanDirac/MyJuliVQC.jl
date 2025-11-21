using Test
using LinearAlgebra
using MyJuliVQC

# Convenient aliases
arity  = MyJuliVQC.arity
kraus  = MyJuliVQC.kraus
qubits = MyJuliVQC.qubits

# ---------------------------
# Utilities
# ---------------------------
I2(T=ComplexF64) = Matrix{T}(I, 2, 2)
X(T=ComplexF64)  = T[0 1; 1 0]
Y(T=ComplexF64)  = T[0 -im; im 0]
Z(T=ComplexF64)  = T[1 0; 0 -1]

@testset "QuantumMap: construction & invariants" begin
    # --- Single-qubit: identity channel (Kraus = {I}) ---
    K = I2()
    M = MyJuliVQC.QuantumMap(1, (K,))  # pos::Int + Tuple{Matrix}
    @test qubits(M) == [1]
    @test arity(M) == 1
    @test length(kraus(M)) == 1
    @test size(kraus(M)[1]) == (2, 2)
    @test eltype(kraus(M)[1]) == ComplexF64

    # --- Single-qubit amplitude damping Kraus (trace preserving) ---
    γ = 0.2
    K0 = ComplexF64[1 0; 0 sqrt(1-γ)]
    K1 = ComplexF64[0 sqrt(γ); 0 0]
    M2 = MyJuliVQC.QuantumMap([1], (K0, K1))  # positions::Vector
    @test qubits(M2) == [1]
    @test arity(M2) == 1
    S = K0' * K0 + K1' * K1
    @test isapprox(S, I2(); atol=1e-10, rtol=0)

    # --- Using Tuple positions overload ---
    M3 = MyJuliVQC.QuantumMap((1,), (K0, K1))
    @test qubits(M3) == [1]
    @test arity(M3) == 1

    # --- Force element type: Float32/ComplexF32 ---
    Kf = Float32[1 0; 0 1]
    M4 = MyJuliVQC.QuantumMap(ComplexF32, 1, (Kf,))   # elty=ComplexF32
    @test eltype(kraus(M4)[1]) == ComplexF32

    # --- Two-qubit identity channel (4×4) ---
    I4 = Matrix{ComplexF64}(I, 4, 4)
    M5 = MyJuliVQC.QuantumMap([1,2], (I4,))  # 4×4 channel acting on two qubits
    @test qubits(M5) == [1, 2]
    @test arity(M5) == 2
    @test size(kraus(M5)[1]) == (4, 4)
end

@testset "QuantumMap: errors & checks" begin
    # --- empty positions ---
    @test_throws ArgumentError MyJuliVQC.QuantumMap(Int[], (I2(),))

    # --- invalid qubit index (<1) ---
    @test_throws ArgumentError MyJuliVQC.QuantumMap([0], (I2(),))

    # --- duplicate qubit indices ---
    @test_throws ArgumentError MyJuliVQC.QuantumMap([1, 1], (I2(),))

    # --- empty Kraus list ---
    @test_throws ArgumentError MyJuliVQC.QuantumMap(1, ())

    # --- non-square Kraus matrix ---
    badK = ComplexF64[1 0 0; 0 1 0]  # 2×3
    @test_throws ArgumentError MyJuliVQC.QuantumMap(1, (badK,))

    # --- dimension mismatch: positions is two qubits but Kraus is 2×2 ---
    @test_throws ArgumentError MyJuliVQC.QuantumMap([1,2], (I2(),))

    # --- trace-preserving check fails: both Kraus = I (sum K†K = 2I) ---
    @test_throws AssertionError MyJuliVQC.QuantumMap(1, (I2(), I2()))

    # --- disable check_tp: allow non-trace-preserving ---
    Mnp = MyJuliVQC.QuantumMap(1, (I2(), I2()); check_tp=false)
    @test qubits(Mnp) == [1]
    @test length(kraus(Mnp)) == 2
end

@testset "QuantumMap: element type normalization" begin
    # input real matrix, expect normalization to ComplexF64
    Kreal = [1.0 0.0; 0.0 1.0]  # Float64
    M = MyJuliVQC.QuantumMap(1, (Kreal,))
    @test eltype(kraus(M)[1]) == ComplexF64

    # mixed element types → unify to ComplexF32, and remain trace-preserving
    α = 1/sqrt(2)
    K1 = α * Float32[1 0; 0 1]            # Float32
    K2 = α * ComplexF64[0 1; 1 0]         # ComplexF64
    M2 = MyJuliVQC.QuantumMap(ComplexF32, 1, (K1, K2))
    @test eltype(kraus(M2)[1]) == ComplexF32
    @test eltype(kraus(M2)[2]) == ComplexF32
    # compute S using internal Kraus operators (ensures consistency)
    Ks = kraus(M2)
    S  = Ks[1]'*Ks[1] + Ks[2]'*Ks[2]
    I2f = Matrix{ComplexF32}(I, 2, 2)

    # Float32 has limited precision → relaxed atol
    @test isapprox(S, I2f; atol=1e-5, rtol=0)
end

@testset "QuantumMap: simple TP channels sanity" begin
    # Phase-damping-like example (coherence damping only)
    γ = 0.35
    K0 = ComplexF64[1 0; 0 sqrt(1-γ)]
    K1 = ComplexF64[0 0; 0 sqrt(γ)]
    M = MyJuliVQC.QuantumMap(1, (K0, K1))
    S = K0' * K0 + K1' * K1
    @test isapprox(S, I2(); atol=1e-10, rtol=0)

    # Depolarizing（common parametrization:K0=√(1-3p/4)I，Ki=√(p/4){X,Y,Z}）
    p = 0.4
    K0 = sqrt(1 - 3p/4) * I2()
    s  = sqrt(p/4)
    K1, K2, K3 = s*X(), s*Y(), s*Z()
    Mdep = MyJuliVQC.QuantumMap(1, (K0, K1, K2, K3))
    Sdep = K0'K0 + K1'K1 + K2'K2 + K3'K3
    @test isapprox(Sdep, I2(); atol=1e-10, rtol=0)
end
