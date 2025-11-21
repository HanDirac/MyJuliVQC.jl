using Test
using LinearAlgebra
using MyJuliVQC

# Assert that `expr` throws `exctype` AND that the error message contains substring `msgsub`
function _throws_with_message(exctype::Type, msgsub::AbstractString, expr::Function)
    ok = false
    try
        expr()
    catch e
        @test e isa exctype
        io = IOBuffer()
        showerror(io, e)
        s = String(take!(io))
        @test occursin(msgsub, s)
        ok = true
    end
    @test ok  # ensure an exception was actually thrown
end
_throws_with_message(expr::Function, exctype::Type, msgsub::AbstractString) =
    _throws_with_message(exctype, msgsub, expr)

@testset "QuantumGate: constructors & validation" begin
    # some commonly used gate matrices
    X = ComplexF64[0 1; 1 0]                # Pauli-X (2×2 unitary)
    I2 = Matrix{ComplexF64}(I, 2, 2)        # 2×2 identity
    I4 = Matrix{ComplexF64}(I, 4, 4)        # 4×4 identity (two-qubit identity gate)
    NotUnitary2 = ComplexF64[1 1; 0 1]      # non-unitary (upper triangular)

    # ---------------- Basic: Vector positions + validate=true (default) ----------------
    @testset "basic (vector positions, validate=true)" begin
        g = QuantumGate([2], X)  # apply X on qubit 2
        @test eltype(MyJuliVQC.matrix(g)) == ComplexF64
        @test MyJuliVQC.arity(g) == 1
        @test MyJuliVQC.qubits(g) == [2]
        @test MyJuliVQC.matrix(g) == X
    end

    # ---------------- Tuple positions overload ---------------
    @testset "tuple positions overload" begin
        # two-qubit gate: identity (4×4)
        g = QuantumGate((1, 3), I4)
        @test MyJuliVQC.arity(g) == 2
        @test MyJuliVQC.qubits(g) == [1, 3]     # order preserved
        @test MyJuliVQC.matrix(g) == I4
    end

    # ---------------- Forced element-type overload ----------------
    @testset "forced eltype overload" begin
        g = QuantumGate(Float32, (1,), I2)      # converted to Float32 and validated for unitarity
        @test eltype(MyJuliVQC.matrix(g)) == Float32
        @test MyJuliVQC.arity(g) == 1
        @test MyJuliVQC.matrix(g) ≈ Matrix{Float32}(I, 2, 2)
    end

    # ---------------- validate=false allows non-unitary ----------------
    @testset "validate=false allows non-unitary" begin
        g = QuantumGate([1], NotUnitary2; validate=false)
        @test MyJuliVQC.arity(g) == 1
        @test MyJuliVQC.matrix(g) == NotUnitary2
    end

    # ---------------- Error branches ----------------
    @testset "error branches" begin
        # 1) positions-related errors
        @test_throws MethodError QuantumGate([], X)  # [] is Vector{Any} → no matching method → MethodError
        _throws_with_message(ArgumentError, "positions must be non-empty") do
        QuantumGate(Int[], X)              # empty position list of correct type
        end
        _throws_with_message(ArgumentError, "qubit indices must be ≥ 1 (1-based)") do
        QuantumGate([0], X)              # not 1-based
        end
        _throws_with_message(ArgumentError, "qubit indices must be unique") do
        QuantumGate([1, 1], I4)            # duplicate indices
        end

        # 2) matrix shape issues
        _throws_with_message(ArgumentError, "gate matrix must be square") do
        QuantumGate([1], zeros(ComplexF64, 2, 3))
        end
        _throws_with_message(ArgumentError, "gate matrix size must be 2^m") do
        QuantumGate([1], zeros(ComplexF64, 3, 3))
        end

        # 3) arity mismatch (1 qubit but given a 4×4 matrix)
        _throws_with_message(ArgumentError, "does not match positions length") do
        QuantumGate([2], I4)   # only 1 qubit but matrix is 4×4
        end

        # 4) Non-unitary with validate=true (default)
#        _throws_with_message(AssertionError, "gate matrix is not unitary") do
#        QuantumGate([1], NotUnitary2)
#        end
    end
end