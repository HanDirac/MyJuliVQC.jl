using Test
using LinearAlgebra
using MyJuliVQC
using MyJuliVQC.Gates

# small helpers
X(T=ComplexF64) = T[0 1; 1 0]
Y(T=ComplexF64) = T[0 -im; im 0]
Z(T=ComplexF64) = T[1 0; 0 -1]
H(T=ComplexF64) = (one(T)/sqrt(T(2))) * T[1 1; 1 -1]
S(T=ComplexF64) = T[1 0; 0 im]
Tg(T=ComplexF64)= T[1 0; 0 exp(one(T)*im* (π/4))]
sqrtX(T=ComplexF64) = (one(T)/2) * T[1+im 1-im; 1-im 1+im]
sqrtY(T=ComplexF64) = (one(T)/2) * T[1+im -1-im; 1+im 1+im]

# helpers to extract gate matrix / qubit positions (names don't conflict with package)
m(g) = MyJuliVQC.matrix(g)
q(g) = MyJuliVQC.qubits(g)

@testset "Gates: 1-qubit standard" begin
    @testset "X/Y/Z/H/S/T" begin
        g = XGate(1)
        @test q(g) == [1]
        @test m(g) == X()

        g = YGate(2)
        @test q(g) == [2]
        @test m(g) == Y()

        g = ZGate(3)
        @test q(g) == [3]
        @test m(g) == Z()

        g = HGate(1)
        @test q(g) == [1]
        @test m(g) ≈ H()

        g = SGate(1)
        @test q(g) == [1]
        @test m(g) == S()

        g = TGate(1)
        @test q(g) == [1]
        @test m(g) == Tg()

        g = sqrtXGate(4)
        @test q(g) == [4]
        @test m(g) ≈ sqrtX()

        g = sqrtYGate(5)
        @test q(g) == [5]
        @test m(g) ≈ sqrtY()

                # --- tuple-position overloads: G((i,)) ---
        g = XGate((1,))
        @test q(g) == [1]
        @test m(g) == X()

        g = YGate((2,))
        @test q(g) == [2]
        @test m(g) == Y()

        g = ZGate((3,))
        @test q(g) == [3]
        @test m(g) == Z()

        g = HGate((1,))
        @test q(g) == [1]
        @test m(g) ≈ H()

        g = SGate((1,))
        @test q(g) == [1]
        @test m(g) == S()

        g = TGate((1,))
        @test q(g) == [1]
        @test m(g) == Tg()

        g = sqrtXGate((4,))
        @test q(g) == [4]
        @test m(g) ≈ sqrtX()

        g = sqrtYGate((5,))
        @test q(g) == [5]
        @test m(g) ≈ sqrtY()
    end
end

@testset "Gates: 1-qubit rotations (matrix & parameterization)" begin
    # non-parametric RxGate/RyGate/RzGate
    θ = π/3
    g = RxGate(1, θ)
    @test q(g) == [1]
    @test m(g) ≈ ComplexF64[
        cos(θ/2)       -im*sin(θ/2);
        -im*sin(θ/2)    cos(θ/2)
    ]

    g = RyGate(2, θ)
    @test q(g) == [2]
    @test m(g) ≈ ComplexF64[
        cos(θ/2)   -sin(θ/2);
        sin(θ/2)    cos(θ/2)
    ]

    g = RzGate(3, θ)
    @test q(g) == [3]
    @test m(g) ≈ ComplexF64[
        exp(-im*θ/2)  0;
        0             exp( im*θ/2)
    ]

        # --- tuple-position overloads: Rx/Ry/Rz((i,), θ) ---
    g = RxGate((1,), θ)
    @test q(g) == [1]
    @test m(g) ≈ ComplexF64[
        cos(θ/2)       -im*sin(θ/2);
        -im*sin(θ/2)    cos(θ/2)
    ]

    g = RyGate((2,), θ)
    @test q(g) == [2]
    @test m(g) ≈ ComplexF64[
        cos(θ/2)   -sin(θ/2);
        sin(θ/2)    cos(θ/2)
    ]

    g = RzGate((3,), θ)
    @test q(g) == [3]
    @test m(g) ≈ ComplexF64[
        exp(-im*θ/2)  0;
        0             exp( im*θ/2)
    ]

    # parameterized RxGate/RyGate/RzGate (isparas = true)
    desc = RxGate(4, θ; isparas=true)
    @test desc isa MyJuliVQC.Gates.ParamOp
    @test length(desc.params) == 1
    @test desc.params[1] == θ
    @test desc.mask == BitVector([true])  
    g2 = desc([π/2])
    @test q(g2) == [4]
    @test m(g2) ≈ ComplexF64[
        cos(π/4)     -im*sin(π/4);
        -im*sin(π/4)  cos(π/4)
    ]

    # RyGate with parameterization flag controlling whether it is actually parametric
    desc2 = RyGate(1, θ; isparas=false)  # not parameterized: returns a concrete QuantumGate
    @test desc2 isa MyJuliVQC.QuantumGate

    desc3 = RzGate(2, θ; isparas=[false])  # fixed parameter: returns a fixed QuantumGate
    @test desc3 isa MyJuliVQC.QuantumGate
    @test q(desc3) == [2]
    @test m(desc3) ≈ ComplexF64[
        exp(-im*θ/2)  0;
        0             exp( im*θ/2)
    ]

        # --- tuple-position overloads: parametric Rx/Ry/Rz((i,), θ; isparas=...) ---
    desc_t = RxGate((4,), θ; isparas=true)
    @test desc_t isa MyJuliVQC.Gates.ParamOp
    @test desc_t.params[1] == θ
    @test desc_t.mask == BitVector([true])
    g_t = desc_t([π/2])
    @test q(g_t) == [4]

    # fixed parameter via isparas=[false] should return concrete QuantumGate
    g_fixed = RzGate((2,), θ; isparas=[false])
    @test g_fixed isa MyJuliVQC.QuantumGate
    @test q(g_fixed) == [2]
    @test m(g_fixed) ≈ ComplexF64[
        exp(-im*θ/2)  0;
        0             exp( im*θ/2)
    ]
end

@testset "Gates: controlled rotations (matrix & parameterization)" begin
    θ = π/5

    # non-parametric CRxGate
    g = CRxGate(1, 2, θ)
    @test q(g) == [1, 2]
    @test m(g) ≈ ComplexF64[
            1 0 0 0;
            0 cos(θ/2) 0 -im*sin(θ/2);
            0 0 1 0;  
            0 -im*sin(θ/2) 0  cos(θ/2)
        ]

        # --- tuple-position overloads: CRx/CRy/CRz((i,j), θ) ---
    g = CRxGate((1, 2), θ)
    @test q(g) == [1, 2]
    @test m(g) ≈ ComplexF64[
        1 0 0 0;
        0 cos(θ/2) 0 -im*sin(θ/2);
        0 0 1 0;
        0 -im*sin(θ/2) 0 cos(θ/2)
    ]

    g = CRyGate((2, 3), θ)
    @test q(g) == [2, 3]
    @test m(g) ≈ ComplexF64[
        1 0 0 0;
        0 cos(θ/2) 0 -sin(θ/2);
        0 0 1 0;
        0 sin(θ/2) 0 cos(θ/2)
    ]

    g = CRzGate((3, 1), θ)
    @test q(g) == [3, 1]
    @test m(g) ≈ ComplexF64[
        1 0 0 0;
        0 exp(-im*θ/2) 0 0;
        0 0 1 0;
        0 0 0 exp(im*θ/2)
    ]
        
    # parameterized CRyGate
    desc = CRyGate(2, 3, θ; isparas=true)
    #@test haskey(desc, :params) && haskey(desc, :mask) && haskey(desc, :build)
    g2 = desc.build([π/3])
    @test q(g2) == [2, 3]
    @test m(g2) ≈ ComplexF64[
        1 0 0 0;
        0 cos(π/6) 0 -sin(π/6);
        0 0 1 0  ;
        0 sin(π/6) 0  cos(π/6)
    ]

        # --- tuple-position overload: parametric CRyGate((i,j), θ; isparas=true) ---
    desc_t = CRyGate((2, 3), θ; isparas=true)
    @test desc_t isa MyJuliVQC.Gates.ParamOp
    @test desc_t.params[1] == θ
    @test desc_t.mask == BitVector([true])
    g_t = desc_t.build([π/3])
    @test q(g_t) == [2, 3]

    # CRzGate with fixed parameter (mask=false)
    desc2 = CRzGate(3, 1, θ; isparas=[false])
    @test desc2 isa MyJuliVQC.QuantumGate
    @test q(desc2) == [3, 1]
    @test m(desc2) ≈ ComplexF64[
        1 0 0 0;
        0 exp(-im*θ/2) 0 0;
        0 0  1  0;
        0 0  0              exp( im*θ/2)
    ]
end

@testset "Gates: 2-qubit fixed" begin
    g = SWAPGate(1, 3)
    @test q(g) == [1, 3]
    @test m(g) == ComplexF64[
        1 0 0 0;
        0 0 1 0;
        0 1 0 0;
        0 0 0 1
    ]

    g = iSWAPGate(2, 4)
    @test q(g) == [2, 4]
    @test m(g) == ComplexF64[
        1 0 0 0;
        0 0 im 0;
        0 im 0 0;
        0 0 0 1
    ]

    g = CNOTGate(1, 2)
    @test q(g) == [1, 2]
    @test m(g) == ComplexF64[
            1 0 0 0;
            0 0 0 1;
            0 0 1 0;
            0 1 0 0]

    g = CZGate(2, 1)
    @test q(g) == [2, 1]
    @test m(g) == ComplexF64[
        1 0 0 0;
        0 1 0 0;
        0 0 1 0;
        0 0 0 -1
    ]

        # --- tuple-position overloads: 2-qubit fixed gates ---
    g = SWAPGate((1, 3))
    @test q(g) == [1, 3]
    @test m(g) == ComplexF64[
        1 0 0 0;
        0 0 1 0;
        0 1 0 0;
        0 0 0 1
    ]

    g = iSWAPGate((2, 4))
    @test q(g) == [2, 4]
    @test m(g) == ComplexF64[
        1 0 0 0;
        0 0 im 0;
        0 im 0 0;
        0 0 0 1
    ]

    g = CNOTGate((1, 2))
    @test q(g) == [1, 2]
    @test m(g) == ComplexF64[
        1 0 0 0;
        0 0 0 1;
        0 0 1 0;
        0 1 0 0
    ]

    g = CZGate((2, 1))
    @test q(g) == [2, 1]
    @test m(g) == ComplexF64[
        1 0 0 0;
        0 1 0 0;
        0 0 1 0;
        0 0 0 -1
    ]
end

@testset "Gates: 3-qubit fixed" begin
    # Toffoli = controlled-controlled-NOT
    g = TOFFOLIGate(1, 2, 3)
    @test q(g) == [1, 2, 3]
    @test m(g) == ComplexF64[
        1 0 0 0 0 0 0 0;
        0 1 0 0 0 0 0 0;
        0 0 1 0 0 0 0 0;
        0 0 0 0 0 0 0 1;
        0 0 0 0 1 0 0 0;
        0 0 0 0 0 1 0 0;
        0 0 0 0 0 0 1 0;  
        0 0 0 1 0 0 0 0  
    ]

   # Fredkin (controlled-SWAP): in the standard definition, when control = 1,
    # the two target qubits are swapped.
    # Expected matrix (basis ordering (c,t1,t2) = 000,001,010,011,100,101,110,111)
F_expected = ComplexF64[
    1 0 0 0 0 0 0 0;
    0 1 0 0 0 0 0 0;
    0 0 1 0 0 0 0 0;
    0 0 0 0 0 1 0 0;
    0 0 0 0 1 0 0 0;
    0 0 0 1 0 0 0 0;
    0 0 0 0 0 0 1 0;  
    0 0 0 0 0 0 0 1 
]

gF = FREDKINGate(1, 2, 3)
@test q(gF) == [1, 2, 3]
@test m(gF) == F_expected

        # --- tuple-position overloads: 3-qubit fixed gates ---
    g = TOFFOLIGate((1, 2, 3))
    @test q(g) == [1, 2, 3]
    @test m(g) == ComplexF64[
        1 0 0 0 0 0 0 0;
        0 1 0 0 0 0 0 0;
        0 0 1 0 0 0 0 0;
        0 0 0 0 0 0 0 1;
        0 0 0 0 1 0 0 0;
        0 0 0 0 0 1 0 0;
        0 0 0 0 0 0 1 0;
        0 0 0 1 0 0 0 0
    ]

    gF2 = FREDKINGate((1, 2, 3))
    @test q(gF2) == [1, 2, 3]
    @test m(gF2) == F_expected

end

@testset "Gates: generic CONTROL/CONTROLCONTROL helpers" begin
    # Use generic helpers to embed Y/X into controlled gates
    g = CONTROLGate(1, 2; U=Y())
    @test q(g) == [1, 2]
    @test m(g) == ComplexF64[
        1 0 0 0;
        0 0 0 -im;
        0 0 1 0;
        0 im 0  0
    ]

    g2 = CONTROLCONTROLGate(1, 2, 3; U=X())
    @test q(g2) == [1, 2, 3]
    @test m(g2) == ComplexF64[
        1 0 0 0 0 0 0 0;
        0 1 0 0 0 0 0 0;
        0 0 1 0 0 0 0 0;
        0 0 0 0 0 0 0 1;
        0 0 0 0 1 0 0 0;
        0 0 0 0 0 1 0 0;
        0 0 0 0 0 0 1 0;  
        0 0 0 1 0 0 0 0 
    ]

        # --- tuple-position overloads for generic helpers ---
    g = CONTROLGate((1, 2); U=Y())
    @test q(g) == [1, 2]
    @test m(g) == ComplexF64[
        1 0 0 0;
        0 0 0 -im;
        0 0 1 0;
        0 im 0  0
    ]

    g2 = CONTROLCONTROLGate((1, 2, 3); U=X())
    @test q(g2) == [1, 2, 3]
    @test m(g2) == ComplexF64[
        1 0 0 0 0 0 0 0;
        0 1 0 0 0 0 0 0;
        0 0 1 0 0 0 0 0;
        0 0 0 0 0 0 0 1;
        0 0 0 0 1 0 0 0;
        0 0 0 0 0 1 0 0;
        0 0 0 0 0 0 1 0;
        0 0 0 1 0 0 0 0
    ]

end

@testset "Gates: FSIM matrix & parameterization" begin
    θ = π/6; ϕ = π/8
    g = FSIMGate(1, 2, θ, ϕ)
    @test q(g) == [1, 2]
    a = cos(θ)
    b = -im * sin(θ)
    c = exp(-im*ϕ)
    @test m(g) ≈ ComplexF64[
        1 0 0 0;
        0 a b 0;
        0 b a 0;
        0 0 0 c
    ]

        # --- tuple-position overloads: FSIMGate((i,j), θ, ϕ) ---
    g = FSIMGate((1, 2), θ, ϕ)
    @test q(g) == [1, 2]
    @test m(g) ≈ ComplexF64[
        1 0 0 0;
        0 a b 0;
        0 b a 0;
        0 0 0 c
    ]

    # parameterized: both parameters are active
    desc = FSIMGate(2, 3, θ, ϕ; isparas=[true, true])
    #@test haskey(desc, :params) && haskey(desc, :mask) && haskey(desc, :build)
    @test desc.params ≈ [θ, ϕ]
    @test desc.mask == BitVector([true, true])
    g2 = desc.build([π/3, π/10])
    @test q(g2) == [2, 3]
    @test m(g2) ≈ ComplexF64[
        1 0 0 0;
        0 cos(π/3) -im*sin(π/3) 0;
        0 -im*sin(π/3) cos(π/3) 0;
        0 0 0 exp(-im*π/10)
    ]

        # --- tuple-position overload: parametric FSIMGate((i,j), ...; isparas=...) ---
    desc_t = FSIMGate((2, 3), θ, ϕ; isparas=[true, true])
    @test desc_t isa MyJuliVQC.Gates.ParamOp
    @test desc_t.params ≈ [θ, ϕ]
    @test desc_t.mask == BitVector([true, true])
    g_t = desc_t.build([π/3, π/10])
    @test q(g_t) == [2, 3]

    # only θ is active
    desc2 = FSIMGate(3, 4, θ, ϕ; isparas=[true, false])
    @test desc2.params ≈ [θ, ϕ]
    @test desc2.mask == BitVector([true, false])
    g3 = desc2.build([π/7, ϕ])  # second parameter stays fixed
    @test q(g3) == [3, 4]
    @test m(g3) ≈ ComplexF64[
        1 0 0 0;
        0 cos(π/7) -im*sin(π/7) 0;
        0 -im*sin(π/7) cos(π/7) 0;
        0 0 0 exp(-im*ϕ)
    ]

    #g4 = desc2.build([π/7, 99999])   # mess up the second parameter
    #@test m(g4) ≈ m(g3)              # result should still be identical to g3
end
