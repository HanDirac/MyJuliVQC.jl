using Test
using LinearAlgebra
using Random
using MyJuliVQC

# Internal helpers needed: packing/writing circuit parameters
# (consistent with src/gradient.jl)
#using MyJuliVQC: flattened_ops

const _kron_many = MyJuliVQC._kron_many

# Heisenberg-1D (constructed with existing QubitsTerm/QubitsOperator)
function heisenberg_1d(L; hz=0.3, J=0.7)
    terms = QubitsTerm[]
    for i in 1:L
        push!(terms, QubitsTerm(i=>"Z"; coeff=hz))
    end
    for i in 1:(L-1)
        push!(terms, QubitsTerm(i=>"X", i+1=>"X"; coeff=J))
        push!(terms, QubitsTerm(i=>"Y", i+1=>"Y"; coeff=J))
        push!(terms, QubitsTerm(i=>"Z", i+1=>"Z"; coeff=J))
    end
    QubitsOperator(terms)
end

# Build a circuit using only the gate constructors provided by the package
function build_circuit(; L=3, depth=4, rng::AbstractRNG=Random.default_rng())
    ops = CircuitElement[]
    for _ in 1:depth
        for i in 1:(L-1)
            push!(ops, CNOTGate(i, i+1))     # built-in 2-qubit gate
        end
        for i in 1:L
            push!(ops, RyGate(i, randn(rng), isparas=true))  # parameterized gate
            push!(ops, RxGate(i, randn(rng), isparas=true))
        end
    end
    QCircuit(ops)
end

# Loss and gradient utilities (directly exported by the package)
const LossExpectationRealSV = MyJuliVQC.LossExpectationRealSV
const LossExpectationRealDM = MyJuliVQC.LossExpectationRealDM
const gradientMJVQC         = MyJuliVQC.gradient
#const apply_ad              = MyJuliVQC.apply_ad

@testset "gradient (SV & DM)" begin
    L  = 3
    sv = StateVector(L)     # |000⟩
    dm = DensityMatrix(L)   # |000⟩⟨000|
    op = heisenberg_1d(L)

    rng = MersenneTwister(7)  
    circ = build_circuit(; L=L, depth=4, rng=rng)

    # Four loss variants
#    lossa(c::QCircuit) = real(expectation(op, apply_ad(c, sv)))
    lossc = LossExpectationRealSV(op, sv)

#    lossd(c::QCircuit) = real(expectation(op, apply_ad(c, dm)))
    lossf = LossExpectationRealDM(op, dm)

#    g_a = gradientMJVQC(lossa, circ)
    g_c = gradientMJVQC(lossc, circ)
#    g_d = gradientMJVQC(lossd, circ)
    g_f = gradientMJVQC(lossf, circ) 

    @test isapprox(g_f, g_c; rtol=1e-5, atol=5e-6) 

    # Verify which gradient method is dispatched (multiple dispatch check)
    let
#        m_gen_sv = @which MyJuliVQC.gradient(lossa, circ)
#        m_gen_dm = @which MyJuliVQC.gradient(lossd, circ)
        m_spc_c  = @which MyJuliVQC.gradient(lossc,  circ)
        m_spc_f  = @which MyJuliVQC.gradient(lossf,  circ)

        # General dispatch should hit gradient(::Function, ::QCircuit)
#        @test occursin("Function", string(m_gen_sv.sig)) && occursin("QCircuit", string(m_gen_sv.sig))
#        @test occursin("Function", string(m_gen_dm.sig)) && occursin("QCircuit", string(m_gen_dm.sig))

        # Specialized dispatch should hit the overloads for state-vector / density-matrix
        @test occursin("LossExpectationRealSV",  string(m_spc_c.sig))
        @test occursin("LossExpectationRealDM",  string(m_spc_f.sig))
    end
end

@testset "additional test for gradient (SV & DM)" begin
    circ = QCircuit([RxGate(1, pi/6, isparas=true), RyGate(2, pi/3, isparas=true)])
    state = StateVector(2)  # |00⟩
    op = QubitsOperator([QubitsTerm(1=>"Z", 2=>"Z"; coeff=1.0)])
    loss_sv = LossExpectationRealSV(op, state)
    grad = gradientMJVQC(loss_sv, circ)
    @test isapprox(grad, [-0.25, -0.75]; rtol=1e-5, atol=5e-6) 
end


