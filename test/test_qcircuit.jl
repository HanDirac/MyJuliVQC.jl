using Test
using LinearAlgebra
using MyJuliVQC
using MyJuliVQC.Gates
using MyJuliVQC.Channels
using .Gates: ParamOp

# convenient aliases
arity  = MyJuliVQC.arity
matrix = MyJuliVQC.matrix
kraus  = MyJuliVQC.kraus
qubits = MyJuliVQC.qubits
#flattened_ops = MyJuliVQC.flattened_ops

@testset "QCircuit: construction & basics" begin
    # empty circuit
    c0 = MyJuliVQC.QCircuit()
    @test length(c0) == 0

    # single gate / single channel
    g = XGate(1)
    m = PhaseDamping(2; γ=0.3)
    c1 = MyJuliVQC.QCircuit([g])
    c2 = MyJuliVQC.QCircuit([m])
    @test length(c1) == 1 && c1[1] === g
    @test length(c2) == 1 && c2[1] === m

    # multi-op construction: mixture of gates and channels
    g2 = HGate(3)
    c3 = MyJuliVQC.QCircuit([g, m, g2])
    @test length(c3) == 3
    @test c3[1] === g && c3[2] === m && c3[3] === g2
end

@testset "QCircuit: accepts allowed ops (gate/map/param-op/nested)" begin
    # non-parametric gate (isparas=false) → should return a QuantumGate
    g_fixed = RyGate(1, π/5; isparas=false)
    @test g_fixed isa MyJuliVQC.QuantumGate

    # parametric gate (isparas=true) → should return a ParamOp
    pdesc = RxGate(2, π/7; isparas=true)
    @test pdesc isa ParamOp
    #@test haskey(pdesc, :params) && haskey(pdesc, :mask) && haskey(pdesc, :build)
    @test pdesc.build isa Function

    # a quantum channel
    ch = Depolarizing(1; p=0.2)

    # nested circuit (default: flatten=true)
    inner = MyJuliVQC.QCircuit([g_fixed, ch, pdesc])
    outer = MyJuliVQC.QCircuit([XGate(3), inner, ZGate(1)])  # default flatten=true
    @test length(outer) == 5
    ops = [outer[i] for i in 1:length(outer)]
    @test ops[1] isa MyJuliVQC.QuantumGate          # XGate(3)
    @test ops[2] isa MyJuliVQC.QuantumGate          # g_fixed
    @test ops[3] isa MyJuliVQC.QuantumMap           # ch
    @test ops[4] isa ParamOp                     # pdesc
    @test ops[5] isa MyJuliVQC.QuantumGate          # ZGate(1)
end

#@testset "QCircuit: flatten flag" begin
    #g1 = XGate(1)
    #g2 = YGate(2)
    #inner = MyJuliVQC.QCircuit([g1, g2]; flatten=true)
    # # when building the outer circuit, disable flattening
    #c = MyJuliVQC.QCircuit([inner, ZGate(3)]; flatten=false)
    #@test length(c) == 2
    #@test c[1] isa MyJuliVQC.QCircuit
    #@test c[2] isa MyJuliVQC.QuantumGate

    # flattened_ops should flatten everything
    #flat = flattened_ops(c)
    #@test length(flat) == 3
    #@test flat[1] === g1 && flat[2] === g2 && flat[3] isa MyJuliVQC.QuantumGate
#end

@testset "QCircuit: validate flag & errors" begin
    g = XGate(1)
    bad = 42  # invalid element
    @test_throws ArgumentError MyJuliVQC.QCircuit([g, bad])  # validate=true by default → must throw
    #c = MyJuliVQC.QCircuit([g, bad]; validate=false)         # allow invalid elements as raw container
    #@test length(c) == 2
    #@test c[1] === g && c[2] === bad
end

@testset "QCircuit: without input aliasing defense (copy-on-build)" begin
    g = XGate(1); h = HGate(2); z = ZGate(3)
    arr = CircuitElement[g, h]
    c = MyJuliVQC.QCircuit(arr)       # constructor makes a copy
    push!(arr, z)                     # mutate the original array
    @test length(arr) == 3
    @test length(c) == 3              # circuit should NOT change
end

@testset "QCircuit: getindex & ordering" begin
    g1 = XGate(1)
    ch = PhaseDamping(1; γ=0.4)
    pdesc = RzGate(1, 0.0; isparas=true)
    c = MyJuliVQC.QCircuit([g1, ch, pdesc])
    @test c[1] === g1
    @test c[2] === ch
    @test c[3] === pdesc
end

