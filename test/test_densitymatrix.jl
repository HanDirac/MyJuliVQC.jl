using Test
using LinearAlgebra
using MyJuliVQC

@testset "DensityMatrix: constructors & validation" begin
    # Basic constructor: |0⟩⟨0|^{⊗n}
    @testset "init |0⟩⟨0|^{⊗n}" begin
        ρ = DensityMatrix(3)
        @test eltype(ρ) == ComplexF64
        @test length(ρ) == 4^3
        @test MyJuliVQC.nqubits(ρ) == 3
        M = MyJuliVQC.matview(ρ)
        @test size(M) == (2^3, 2^3)
        @test M[1,1] == 1
        @test isapprox(LinearAlgebra.tr(M), 1.0)
    end

    # type parameter
    @testset "type parameter" begin
        ρf = DensityMatrix(Float32, 2)
        @test eltype(ρf) == Float32
        @test MyJuliVQC.nqubits(ρf) == 2
    end

    # n = 0 corner case
    @testset "n=0 corner case" begin
        ρ0 = DensityMatrix(0)
        @test MyJuliVQC.nqubits(ρ0) == 0
        @test length(ρ0) == 1
        @test MyJuliVQC.matview(ρ0)[1,1] == 1
    end

    # Construct from a matrix (already a valid density matrix)
    @testset "from matrix (valid)" begin
        A = zeros(ComplexF64, 4, 4); A[1,1] = 0.5; A[2,2] = 0.5
        ρA = DensityMatrix(A; validate=true)
        @test MyJuliVQC.nqubits(ρA) == 2
        @test isapprox(LinearAlgebra.tr(MyJuliVQC.matview(ρA)), 1.0)
    end

    # Construct from a vector (column-major); must be consistent with matrix form
    @testset "from vectorized (column-major)" begin
        Mcol = reshape(collect(1:16), 4, 4) # column-major: 1 5 9 13; 2 6 10 14; ...
        v = vec(Mcol)
        ρv = DensityMatrix(v)
        @test MyJuliVQC.matview(ρv) == Mcol
    end

    # normalize_trace: scale trace to 1
    @testset "normalize_trace + validate" begin
        B = zeros(ComplexF64, 2, 2); B[1,1] = 2.0
        ρB = DensityMatrix(B; normalize_trace=true, validate=true)
        @test isapprox(LinearAlgebra.tr(MyJuliVQC.matview(ρB)), 1.0)
    end

    # error branches
    @testset "error branches" begin
        @test_throws ArgumentError DensityMatrix(zeros(3,4))     # non-square matrix
        @test_throws ArgumentError DensityMatrix(zeros(3,3))     # dimension not 2^n
        @test_throws AssertionError DensityMatrix([1.0, 0.0, 0.0])# vector length ≠ 2^(2n) (2025-08-19: currently this test does not pass)

        NH = [1 1+0im 0 0;
              0 1     0 0;
              0 0     0 0;
              0 0     0 0]
        @test_throws AssertionError DensityMatrix(NH; validate=true)  # not Hermitian

        NPSD = zeros(ComplexF64, 2, 2); NPSD[1,1] = 1.0; NPSD[2,2] = -0.1
        @test_throws AssertionError DensityMatrix(NPSD; normalize_trace=true, validate=true) # not PSD

        C = zeros(ComplexF64, 2, 2); C[1,1] = 0.3; C[2,2] = 0.3  # trace = 0.6
        @test_throws AssertionError DensityMatrix(C; validate=true)  # trace ≠ 1
        ρC = DensityMatrix(C; normalize_trace=true, validate=true)   # normalize then validate
        @test isapprox(LinearAlgebra.tr(MyJuliVQC.matview(ρC)), 1.0)
    end
end
