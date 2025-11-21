using Test
using LinearAlgebra
using MyJuliVQC
using MyJuliVQC.Gates
using MyJuliVQC.Channels

# For readability
ap = MyJuliVQC.active_parameters

@testset "active_parameters: basics" begin
    # Empty circuit
    c0 = QCircuit()
    ps0 = ap(c0)
    @test ps0 isa Vector
    @test length(ps0) == 0

    # All fixed (non-parameterized) gates
    c1 = QCircuit([XGate(1), ZGate(2), HGate(3)])
    ps1 = ap(c1)
    @test length(ps1) == 0
end

@testset "active_parameters: simple param gates" begin
    θ1 = π/3
    θ2 = π/5

    # RxGate is parameterized; RyGate with isparas=false returns a fixed QuantumGate
    op1 = RxGate(1, θ1; isparas=true)
    op2 = RyGate(2, θ2; isparas=false)

    c = QCircuit([op1, op2])
    ps = ap(c)

    @test length(ps) == 1
    @test isapprox.(ps, [θ1]) |> all
end

@testset "active_parameters: nested circuits & flattening" begin
    θa = π/7
    ϕa = π/11
    # FSIM has two parameters; only the first one is active
    opA = FSIMGate(1, 2, θa, ϕa; isparas=[true, false])
    # Another param gate
    θb = π/4
    opB = RzGate(3, θb; isparas=true)

    # Nested structure: top = [XGate, nested], where nested = [opA, [opB]]
    inner = QCircuit([opB])                    # first-level nesting
    nested = QCircuit([opA, inner])            # second-level nesting
    top = QCircuit([XGate(1), nested])         # top level

    ps = ap(top)

    # Expected order: FSIM's active θa first, then opB's θb
    @test length(ps) == 2
    @test isapprox.(ps, [θa, θb]) |> all
end

@testset "active_parameters: multi-parameter mask filtering & order" begin
    # FSIM: both parameters are active
    θ = π/9
    ϕ = π/10
    op = FSIMGate(2, 3, θ, ϕ; isparas=[true, true])

    # Another single-parameter gate but mask=false → treated as fixed gate
    op_fixed = RxGate(1, π/6; isparas=false)

    c = QCircuit([op_fixed, op])
    ps = ap(c)

    @test length(ps) == 2
    # Must follow internal parameter order: [θ, ϕ]
    @test isapprox.(ps, [θ, ϕ]) |> all
end

@testset "active_parameters: combination with fixed gates in between" begin
    # param -> fixed -> param
    t1 = π/8
    t2 = π/12
    p1 = RyGate(1, t1; isparas=true)
    fix = CNOTGate(1, 2)
    p2 = RzGate(2, t2; isparas=true)

    c = QCircuit([p1, fix, p2])
    ps = ap(c)

    @test length(ps) == 2
    @test isapprox.(ps, [t1, t2]) |> all
end

@testset "active_parameters: mask with false entries is ignored" begin
    # Single-parameter gate but mask=[false] → becomes a fixed QuantumGate
    θ = π/5
    desc = RzGate(1, θ; isparas=[false])  # This returns a non-parametric QuantumGate
    c = QCircuit([desc])
    ps = ap(c)
    @test length(ps) == 0
end

@testset "active_parameters: channels contribute no params" begin
    # Channels alone have no parameters
    c0 = QCircuit([AmplitudeDamping(1; γ=0.2)])
    @test length(ap(c0)) == 0

    # Mix: param gate + channels + fixed gate
    θ = π/6
    g1 = RxGate(1, θ; isparas=true)           # 1 active parameter
    ch1 = PhaseDamping(1; γ=0.3)          # no parameters
    ch2 = Depolarizing(2; p=0.25)         # no parameters
    g2 = HGate(2)                         # fixed gate

    c1 = QCircuit([g1, ch1, ch2, g2])
    ps1 = ap(c1)
    @test length(ps1) == 1
    @test isapprox.(ps1, [θ]) |> all
end

@testset "active_parameters: nested circuits with channels" begin
    # Inner: channel + param gate
    θ_in  = π/7
    inner = QCircuit([
        AmplitudeDamping(1; γ=0.1),
        FSIMGate(1, 2, θ_in, π/11; isparas=[true, false])   # only θ active
    ])

    # Outer: channel + param gate
    θ_out = π/5
    outer = QCircuit([
        PhaseDamping(2; γ=0.4),
        inner,
        RzGate(3, θ_out; isparas=true)       
    ])

    ps = ap(outer)
    # Flattening order: first FSIM's θ_in, then the outer RzGate's θ_out
    @test length(ps) == 2
    @test isapprox.(ps, [θ_in, θ_out]) |> all
end

@testset "active_parameters: channels-only circuit" begin
    c = QCircuit([
        Depolarizing(1; p=0.1),
        PhaseDamping(2; γ=0.2),
        AmplitudeDamping(3; γ=0.05)
    ])
    @test length(ap(c)) == 0
end
