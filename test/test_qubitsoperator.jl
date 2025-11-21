using Test
using LinearAlgebra
using MyJuliVQC  
using MyJuliVQC: positions, oplist, coeff, simplify!, numeric_eltype

# Helper: capture the string output of show(io, x)
capture_show(x) = (buf = IOBuffer(); show(buf, x); String(take!(buf)))

@testset "QubitsOperator basic API" begin
    # Basic construction: empty, varargs, and Vector
    op_empty = QubitsOperator()
    @test isempty(op_empty)
    @test length(op_empty) == 0

    t1 = QubitsTerm(1=>"X"; coeff=0.5)
    t2 = QubitsTerm(2=>"Z"; coeff=1.0)
    t3 = QubitsTerm(3=>"Y"; coeff=2.0)  # Note: Y has imaginary entries

    op_var = QubitsOperator(t1, t2)
    @test !isempty(op_var)
    @test length(op_var) == 2
    @test op_var[1] === t1
    @test op_var[2] === t2
    # iterator should work
    @test collect(op_var) == [t1, t2]

    op_vec = QubitsOperator([t2, t3])
    @test length(op_vec) == 2
    @test op_vec[1] === t2
    @test op_vec[2] === t3

    # push! behavior and chainable return
    op_push = QubitsOperator()
    push!(op_push, t1)
    @test length(op_push) == 1
    push!(op_push, t2)
    @test length(op_push) == 2
    @test op_push[1] === t1 && op_push[2] === t2

    # copy should be independent
    op_copy = copy(op_var)
    @test op_copy !== op_var
    @test op_copy.terms !== op_var.terms   # ensure internal vector is copied
    @test op_copy[1] === t1 && op_copy[2] === t2

    # show: check only key fragments
    shown = capture_show(op_var)
    @test occursin("QubitsOperator with 2 term(s):", shown)
    @test occursin("positions=", shown)
    @test occursin("coeff=", shown)
end

@testset "QubitsOperator algebra" begin
    t1 = QubitsTerm(1=>"X"; coeff=1.0)
    t2 = QubitsTerm(2=>"Z"; coeff=2.0)
    t3 = QubitsTerm(1=>"Y", 3=>"Z"; coeff=1.0 + 0.5im)  
    op = QubitsOperator(t1, t2, t3)

    # adjoint: X/Z/Y are Hermitian, coefficients conjugated
    op_dag = adjoint(op)
    @test length(op_dag) == 3
    @test positions(op_dag[1]) == positions(t1) == [1]
    @test positions(op_dag[2]) == positions(t2) == [2]
    @test positions(op_dag[3]) == positions(t3) == [1,3]  
    @test coeff(op_dag[1]) == conj(coeff(t1))
    @test coeff(op_dag[2]) == conj(coeff(t2))
    @test coeff(op_dag[3]) == conj(coeff(t3))
    # oplist matrices for X/Z/Y should be identical (Hermitian)
    @test all(oplist(op_dag[1])[1] .== oplist(t1)[1])
    @test all(oplist(op_dag[2])[1] .== oplist(t2)[1])
    @test all(oplist(op_dag[3])[1] .== oplist(t3)[1])
    @test all(oplist(op_dag[3])[2] .== oplist(t3)[2])

    # Scalar multiplication (left and right)
    λ = 3.0
    opL = λ * op
    opR = op * λ
    @test length(opL) == length(op)
    @test length(opR) == length(op)
    @test coeff(opL[1]) == λ * coeff(t1)
    @test coeff(opL[2]) == λ * coeff(t2)
    @test coeff(opR[1]) == λ * coeff(t1)
    @test coeff(opR[2]) == λ * coeff(t2)

    # Addition (concatenation of terms)
    t3 = QubitsTerm(3=>"Y"; coeff=1.5)
    A = QubitsOperator(t1)
    B = QubitsOperator(t2, t3)
    C = A + B
    @test length(C) == 3
    @test C[1] === t1 && C[2] === t2 && C[3] === t3
end

@testset "QubitsOperator simplify!" begin
    # Two identical terms (input order differs) must be merged
    term_a = QubitsTerm([1,2], ["X","Z"], 0.5)
    term_b = QubitsTerm([2,1], ["Z","X"], 1.5) # normalized to same canonical form

    # normalized positions must match
    @test positions(term_a) == [1,2]
    @test positions(term_b) == [1,2]
    # each site's operator matrices must match (X, Z)
    @test all(oplist(term_a)[1] .== oplist(term_b)[1])
    @test all(oplist(term_a)[2] .== oplist(term_b)[2])

    op = QubitsOperator(term_a, term_b)
    @test length(op) == 2

    simplify!(op)
    @test length(op) == 1  
    only_term = op[1]
    @test positions(only_term) == [1,2]
    @test coeff(only_term) == 2.0

    # Using string constructors must also merge properly
    op2 = QubitsOperator(QubitsTerm(1=>"X", 2=>"Z"; coeff=1.0),
                         QubitsTerm(1=>"X", 2=>"Z"; coeff=3.0))
    simplify!(op2)
    @test length(op2) == 1
    @test coeff(op2[1]) == 4.0
end

@testset "QubitsOperator eltype" begin
    @test numeric_eltype(QubitsOperator()) == ComplexF64

    op_real = QubitsOperator(QubitsTerm(1=>"X"; coeff=1.0))
    @test numeric_eltype(op_real) == Float64 || numeric_eltype(op_real) == ComplexF64

    op_complex_by_matrix = QubitsOperator(QubitsTerm(2=>"Y"; coeff=1.0))
    @test numeric_eltype(op_complex_by_matrix) == ComplexF64

    op_complex_by_coeff = QubitsOperator(QubitsTerm(3=>"Z"; coeff=1.0 + 0.1im))
    @test numeric_eltype(op_complex_by_coeff) == ComplexF64
end
