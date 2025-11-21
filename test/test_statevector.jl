using Test
using MyJuliVQC
using LinearAlgebra

@testset "StateVector constructors" begin
    # |0…0⟩
    ψ = StateVector(3)
    @test eltype(ψ) == ComplexF64
    @test length(ψ) == 2^3
    @test ψ[1] == 1
    @test all(ψ.data[2:end] .== 0)
    @test MyJuliVQC.nqubits(ψ) == 3

    # Force element type
    ψr = StateVector(Float64, 2)
    @test eltype(ψr) == Float64
    @test ψr[1] == 1.0

    # Construct from array (unnormalized)
    v = [1.0, 1.0, 0.0, 0.0]
    ψv = StateVector(v)
    @test length(ψv) == 4
    @test ψv[1] == 1.0

    # Construct from array (normalized)
    ψvn = StateVector(v; normalize=true)
    @test isapprox(norm(ψvn.data), 1.0)

    # Type conversion constructor
    ψc = StateVector(ComplexF32, v)
    @test eltype(ψc) == ComplexF32

    # Zero-qubit state (length 1)
    ψ0 = StateVector(0)
    @test length(ψ0) == 1
    @test MyJuliVQC.nqubits(ψ0) == 0

    # Length not a power of 2 should throw an error
    @test_throws AssertionError StateVector([1.0, 0.0, 0.0])
    # Negative n should throw an error
    @test_throws AssertionError StateVector(-1)
end