using Test
using LinearAlgebra
using Random
using MyJuliVQC

# Let qubit 1 be the global LSB (fastest-changing index), and construct
# |ψ⟩ = q4 ⊗ q3 ⊗ q2 ⊗ q1
kron4(q1, q2, q3, q4) = kron(kron(kron(q4, q3), q2), q1)
normed(v) = v ./ norm(v)

@testset "measure!: pure state on qubit 2" begin
    # -- Construct |ψ₁⟩ --
    # For the 2nd qubit, choose the single-qubit state [√3, 1].
    # After normalization, P(1)=1/4 and P(0)=3/4.
    q1 = normed([5.0, 6.0im])                # Other qubits chosen arbitrarily (nonzero)
    q2 = normed([sqrt(3.0), 1.0])            # Key: probability distribution of qubit 2
    q3 = normed([3.0, 4.0im])
    q4 = normed([1.0, 2.0im])

    ψvec = kron4(q1, q2, q3, q4) |> normed
    ψ    = MyJuliVQC.StateVector(ψvec)       # Already normalized

    # Fix the RNG seed for reproducibility
    rng = MersenneTwister(1234)
    out, prob = MyJuliVQC.measure!(ψ, 2; rng=rng)

    # Expected outcome: either (0, 0.75) or (1, 0.25)
    if out == 0
        @test isapprox(prob, 0.75; atol=1e-12, rtol=0)
    else
        @test out == 1
        @test isapprox(prob, 0.25; atol=1e-12, rtol=0)
    end

    # Check that the post-measurement state collapses correctly
    if out == 0
        # Outcome 0: |ψ'⟩ = (I⊗I⊗|0⟩⟨0|⊗I)|ψ⟩ / √P(0)
        P0 = [1.0 0.0; 0.0 0.0]
        ψ_expect = normed(kron4(q1, P0*q2, q3, q4))
        @test isapprox(ψ.data, ψ_expect; atol=1e-12, rtol=0)
    else
        # Outcome 1: |ψ'⟩ = (I⊗I⊗|1⟩⟨1|⊗I)|ψ⟩ / √P(1)
        P1 = [0.0 0.0; 0.0 1.0]
        ψ_expect = normed(kron4(q1, P1*q2, q3, q4))
        @test isapprox(ψ.data, ψ_expect; atol=1e-12, rtol=0)
    end

    # A repeated measurement on the same qubit must yield the same outcome
    # with probability 1 (post-measurement consistency)
    out2, prob2 = MyJuliVQC.measure!(ψ, 2; rng=rng)
    @test out2 == out
    @test isapprox(prob2, 1.0; atol=1e-12, rtol=0)
end

@testset "measure!: density matrix mixture on qubit 2" begin
    # -- Construct |ψ₁⟩ (same as above) --
    q1a = normed([5.0, 6.0im])
    q2a = normed([sqrt(3.0), 1.0])           # P1=1/4
    q3a = normed([3.0, 4.0im])
    q4a = normed([1.0, 2.0im])
    ψ1   = normed(kron4(q1a, q2a, q3a, q4a))

    # -- Construct |ψ₂⟩ --
    # Let qubit 2 be [√2, 1], normalized → P(1)=1/3, P(0)=2/3
    q1b = normed([11.0, 12.0im])
    q2b = normed([sqrt(2.0), 1.0])           # Key: P1 = 1/3
    q3b = normed([9.0, 10.0im])
    q4b = normed([7.0, 8.0im])
    ψ2   = normed(kron4(q1b, q2b, q3b, q4b))

    # ρ = 0.3|ψ₁⟩⟨ψ₁| + 0.7|ψ₂⟩⟨ψ₂|
    ρmat = 0.3*(ψ1*ψ1') + 0.7*(ψ2*ψ2')
    ρ    = MyJuliVQC.DensityMatrix(ρmat)

    # Theoretical measurement probabilities on qubit 2:
    # P(1) = 0.3*(1/4) + 0.7*(1/3) = 37/120 ≈ 0.30833
    # P(0) = 1 - P(1) = 83/120 ≈ 0.69167
    p1_expect = 37/120
    p0_expect = 83/120

    rng = MersenneTwister(5678)
    out, prob = MyJuliVQC.measure!(ρ, 2; rng=rng)

    if out == 0
        @test isapprox(prob, p0_expect; atol=1e-12, rtol=0)
    else
        @test out == 1
        @test isapprox(prob, p1_expect; atol=1e-12, rtol=0)
    end

    # Check post-measurement collapse
    if out == 0
        #  ρ' = (I⊗I⊗|0⟩⟨0|⊗I) ρ (I⊗I⊗|0⟩⟨0|⊗I) / P(0)
        P0 = [1.0 0.0; 0.0 0.0]
        proj = kron4(I(2), P0, I(2), I(2))
        ρ_expect = proj * ρmat * proj' / p0_expect
        @test isapprox(reshape(ρ.data, 16, 16), ρ_expect; atol=1e-12, rtol=0)
    else
        #  ρ' = (I⊗I⊗|1⟩⟨1|⊗I) ρ (I⊗I⊗|1⟩⟨1|⊗I) / P(1)
        P1 = [0.0 0.0; 0.0 1.0]
        proj = kron4(I(2), P1, I(2), I(2))
        ρ_expect = proj * ρmat * proj' / p1_expect
        @test isapprox(reshape(ρ.data, 16, 16), ρ_expect; atol=1e-12, rtol=0)
    end

    # Repeated measurement must now be deterministic
    out2, prob2 = MyJuliVQC.measure!(ρ, 2; rng=rng)
    @test out2 == out
    @test isapprox(prob2, 1.0; atol=1e-12, rtol=0)
end


# ==================== Batch measurements (always from the original state) ====================

# Reconstruct |ψ₁⟩ and ρ (0.3|ψ₁⟩⟨ψ₁| + 0.7|ψ₂⟩⟨ψ₂|)
# Note: must rebuild here because previous @testsets introduce local scope.
build_psai1_statevector() = begin
    normed(v) = v ./ norm(v)
    kron4(q1, q2, q3, q4) = kron(kron(kron(q4, q3), q2), q1)

    q1 = normed([5.0, 6.0im])
    q2 = normed([sqrt(3.0), 1.0])   # qubit 2: P(1)=1/4
    q3 = normed([3.0, 4.0im])
    q4 = normed([1.0, 2.0im])

    ψvec = normed(kron4(q1, q2, q3, q4))
    MyJuliVQC.StateVector(ψvec)
end

build_rho_mixture() = begin
    normed(v) = v ./ norm(v)
    kron4(q1, q2, q3, q4) = kron(kron(kron(q4, q3), q2), q1)

    # |ψ₁⟩
    q1a = normed([5.0, 6.0im])
    q2a = normed([sqrt(3.0), 1.0])     # P1=1/4
    q3a = normed([3.0, 4.0im])
    q4a = normed([1.0, 2.0im])
    ψ1   = normed(kron4(q1a, q2a, q3a, q4a))

    # |ψ₂⟩
    q1b = normed([11.0, 12.0im])
    q2b = normed([sqrt(2.0), 1.0])     # P1=1/3
    q3b = normed([9.0, 10.0im])
    q4b = normed([7.0, 8.0im])
    ψ2   = normed(kron4(q1b, q2b, q3b, q4b))

    ρmat = 0.3*(ψ1*ψ1') + 0.7*(ψ2*ψ2')
    MyJuliVQC.DensityMatrix(ρmat)
end

# ---------- 1) Perform 120 measurements on |ψ₁⟩ (each trial starts from a fresh copy) ----------
begin
    psai_1_orig = build_psai1_statevector()
    meas_qubit  = 2
    rng = MersenneTwister(123456)   

    results = Int[]
    for t in 1:120
        ψ_copy = MyJuliVQC.StateVector(copy(psai_1_orig.data))  
        out, _ = MyJuliVQC.measure!(ψ_copy, meas_qubit; rng=rng)
        push!(results, out)
    end

    println("\n=== 120 independent measurements on psai_1 (always from the original state) ===")
    println("Sequence:")
    println(results)
    println("Counts: #0 = $(count(==(0), results)), #1 = $(count(==(1), results))")
end

# ---------- 2) Perform 120 measurements on ρ (each trial starts from a fresh copy) ----------
begin
    rho_orig  = build_rho_mixture()
    meas_qubit = 2
    rng = MersenneTwister(789012)   

    results = Int[]
    for t in 1:120
        ρ_copy = MyJuliVQC.DensityMatrix(copy(rho_orig.data))   
        out, _ = MyJuliVQC.measure!(ρ_copy, meas_qubit; rng=rng)
        push!(results, out)
    end

    println("\n=== 120 independent measurements on rho (always from the original state) ===")
    println("Sequence:")
    println(results)
    println("Counts: #0 = $(count(==(0), results)), #1 = $(count(==(1), results))")
end
