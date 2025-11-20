using Test
using LinearAlgebra
using MyJuliVQC
using MyJuliVQC.Gates
using MyJuliVQC.Channels

ap = MyJuliVQC.active_parameters
rp! = MyJuliVQC.reset_parameters!

# Small helper: check whether an object is a "parametric operation descriptor"
_is_paramop(x) = x isa MyJuliVQC.Gates.ParamOp

# Check whether two ParamOp objects are equivalent (semantically), not literally equal
function paramops_equivalent(a, b; atol_par = sqrt(eps(Float64)), atol_mat = 1e-10)
    _is_paramop(a) && _is_paramop(b) || return false
    # masks must match
    a.mask == b.mask || return false
    # parameters approximately match (promote to Float64 before comparison to avoid eltype mismatch)
    all(isapprox.(Float64.(a.params), Float64.(b.params); atol=atol_par, rtol=0)) || return false
    # gates built by `build` must match (same qubit positions and matrices)
    ga = a(a.params)  # ParamOp is callable
    gb = b(b.params)
    ga isa MyJuliVQC.QuantumGate && gb isa MyJuliVQC.QuantumGate ||
        return false
    MyJuliVQC.qubits(ga) == MyJuliVQC.qubits(gb) ||
        return false
    isapprox(MyJuliVQC.matrix(ga), MyJuliVQC.matrix(gb); atol=atol_mat, rtol=0)
end

@testset "reset_parameters!: basic single-parameter op" begin
    θ = π/3
    op = RxGate(1, θ; isparas=true)       # NamedTuple
    c  = QCircuit([op])

    # initial parameters
    @test isapprox.(ap(c), [θ]) |> all

    # reset to new value
    θ_new = π/5
    rp!(c, [θ_new])

    # still a "parametric operation" (ParamOp); active_parameters should see the new value
    @test length(c) == 1
    @test _is_paramop(c[1])
    @test isapprox.(ap(c), [θ_new]) |> all
end

@testset "reset_parameters!: multi-parameter with mask" begin
    θ = π/7
    ϕ = π/11
    # only the first parameter is trainable
    op = FSIMGate(1, 2, θ, ϕ; isparas=[true, false])
    c  = QCircuit([op])

    @test isapprox.(ap(c), [θ]) |> all

    θ_new = π/9
    rp!(c, [θ_new])

    # still exactly 1 active parameter, updated to the new value
    @test length(ap(c)) == 1
    @test isapprox.(ap(c), [θ_new]) |> all
    @test _is_paramop(c[1])
end

@testset "reset_parameters!: both params variable, and order preserved" begin
    θ = π/10
    ϕ = π/12
    op = FSIMGate(2, 3, θ, ϕ; isparas=[true, true])
    c  = QCircuit([op])

    @test isapprox.(ap(c), [θ, ϕ]) |> all

    θ2, ϕ2 = π/8, π/6
    rp!(c, [θ2, ϕ2])

    # order must remain [θ, ϕ]
    @test isapprox.(ap(c), [θ2, ϕ2]) |> all
    @test _is_paramop(c[1])
end

@testset "reset_parameters!: interleave with fixed gates and channels" begin
    # trainable + fixed gate + channel + trainable
    t1 = π/8
    t2 = π/12
    p1 = RyGate(1, t1; isparas=true)            # 1 param
    fix = CNOTGate(1, 2)                    # fixed
    ch  = Depolarizing(2; p=0.2)            # channel
    p2 = RzGate(2, t2; isparas=true)            # 1 param

    c = QCircuit([p1, fix, ch, p2])
    @test isapprox.(ap(c), [t1, t2]) |> all

    new_t1 = π/5
    new_t2 = π/7
    rp!(c, [new_t1, new_t2])

    @test isapprox.(ap(c), [new_t1, new_t2]) |> all
    # ensure both parametric positions are still ParamOp
    idxs = [i for i in 1:length(c) if _is_paramop(c[i])]
    @test length(idxs) == 2

    # Check semantic equivalence with the expected ParamOps (instead of literal equality)
    expect1 = RyGate(1, new_t1; isparas=true)
    expect4 = RzGate(2, new_t2; isparas=true)
    @test paramops_equivalent(c[1], expect1)
    @test paramops_equivalent(c[4], expect4)

    # Also confirm that the middle two elements keep their types and positions
    @test c[2] isa MyJuliVQC.QuantumGate
    @test MyJuliVQC.qubits(c[2]) == [1, 2]     # CNOT(1,2)
    @test c[3] isa MyJuliVQC.QuantumMap
    @test MyJuliVQC.qubits(c[3]) == [2]        # Depolarizing(2; p=0.2)

        
end

@testset "reset_parameters!: nested circuits (constructor flattens), still works" begin
    θ_in  = π/7
    inner = QCircuit([FSIMGate(1, 2, θ_in, π/13; isparas=[true, false])])  # 1 param
    θ_out = π/5
    outer = QCircuit([inner, RzGate(3, θ_out; isparas=true)])                  # +1 param
    c     = QCircuit([XGate(1), outer])  # default flatten=true; outer gets flattened

    @test isapprox.(ap(c), [θ_in, θ_out]) |> all

    θ_in2, θ_out2 = π/4, π/9
    rp!(c, [θ_in2, θ_out2])

    @test isapprox.(ap(c), [θ_in2, θ_out2]) |> all
    # should still have exactly 2 parametric operations
    @test sum(_is_paramop, c.operations) == 2

    # Check semantic equivalence to the “ideal” description (after flattening they sit at indices 2 and 3)
    expect_fsim = FSIMGate(1, 2, θ_in2, π/13; isparas=[true, false])
    expect_rz   = RzGate(3, θ_out2; isparas=true)

    @test paramops_equivalent(c[2], expect_fsim)
    @test paramops_equivalent(c[3], expect_rz)

    # Also confirm the first element is still the fixed XGate(1)
    @test c[1] isa MyJuliVQC.QuantumGate
    @test MyJuliVQC.qubits(c[1]) == [1]

end

@testset "reset_parameters!: insufficient / excessive parameters" begin
    # two-parameter gate
    op = FSIMGate(1, 2, π/10, π/12; isparas=[true, true])
    c  = QCircuit([op])

    # too few parameters -> error
    @test_throws ArgumentError rp!(c, [π/3])

    # too many parameters -> error
    @test_throws ArgumentError rp!(c, [π/3, π/4, π/5])

    # exactly 2 parameters -> OK
    @test rp!(c, [π/6, π/7]) === c
    @test isapprox.(ap(c), [π/6, π/7]) |> all
end

@testset "reset_parameters!: multiple calls are allowed (remain param-ops)" begin
    p1 = RxGate(1, π/9; isparas=true)      # 1 param
    p2 = RyGate(2, π/8; isparas=true)      # 1 param
    c  = QCircuit([p1, p2])

    @test isapprox.(ap(c), [π/9, π/8]) |> all
    rp!(c, [π/4, π/6])
    @test isapprox.(ap(c), [π/4, π/6]) |> all

    # call again
    rp!(c, [π/3, π/5])
    @test isapprox.(ap(c), [π/3, π/5]) |> all

    # still parametric operations
    @test sum(_is_paramop, c.operations) == 2
end

@testset "reset_parameters!: type conversion follows original params eltype" begin
    # Build a gate whose params are stored as ComplexF32 (via T=ComplexF32)
    # Use RyGate via the general param-op path
    θ = Float32(π/6)
    op = RyGate(1, θ; T=ComplexF32, isparas=true)   # internal params should be of ComplexF32-based type
    c  = QCircuit([op])

    # Provide new parameters as Float64; reset_parameters! should convert to eltype(p0)
    θ_new = π/7   # Float64
    rp!(c, [θ_new])

    ps = ap(c)
    @test length(ps) == 1
    # allow small numerical error
    @test isapprox(ps[1], θ_new; atol=1e-6)
end
