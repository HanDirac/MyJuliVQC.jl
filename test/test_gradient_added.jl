using Test
using MyJuliVQC

const gradientMJVQC = MyJuliVQC.gradient

#log_path = "grad_desc_log_noisy.txt"
#io = open(log_path, "w")   

#function logprintln(io::IO, args...)
#    println(args...)           # stdout
#    println(io, args...)       # file
#    flush(io)                  
#end
@testset "Compare noisy gradient with PennyLane" begin

# -----------------------------
# 1. Problem setup
# -----------------------------
L      = 3                    # number of qubits
depth  = 2                    # circuit depth
ψ0     = DensityMatrix(L)       # initial state |000⟩⟨000|

# Simple 1D Heisenberg Hamiltonian as QubitsOperator:
#   H = ∑_i hz * Z_i + ∑_i J * (X_i X_{i+1} + Y_i Y_{i+1} + Z_i Z_{i+1})
function heisenberg_1d(L; hz = 1.0, J = 1.0)
    terms = QubitsTerm[]

    # Local Z fields
    for i in 1:L
        push!(terms, QubitsTerm(i => "Z"; coeff = hz))
    end

    # Nearest-neighbor interactions: X X + Y Y + Z Z
    for i in 1:L-1
        push!(terms, QubitsTerm(i => "X", i+1 => "X"; coeff = J))
        push!(terms, QubitsTerm(i => "Y", i+1 => "Y"; coeff = J))
        push!(terms, QubitsTerm(i => "Z", i+1 => "Z"; coeff = J))
    end

    return QubitsOperator(terms)
end

H = heisenberg_1d(L)     # user-defined Hamiltonian as QubitsOperator

# -----------------------------
# 2. Build a variational circuit
# -----------------------------
circuit = QCircuit()
noise_p = 0.1   # depolarizing strength
noise_gamma = 0.3   # amplitude damping strength
noise_gamma_phase = 0.5   # phase damping strength

# First layer: local rotations
for i in 1:L
    push!(circuit, RzGate(i, pi/7; isparas = true))
    push!(circuit, RyGate(i, pi/7; isparas = true))
    push!(circuit, RzGate(i, pi/7; isparas = true))
end

# Repeated entangling + rotations
for l in 1:depth
    # entangling chain
    for i in 1:L-1
        push!(circuit, CNOTGate(i, i+1))
    end
    # another layer of local rotations
    for i in 1:L
        push!(circuit, RzGate(i, pi/7; isparas = true))
        push!(circuit, RxGate(i, pi/7; isparas = true))
        push!(circuit, RzGate(i, pi/7; isparas = true))
    end
    # add noise
    for i in 1:L
        push!(circuit, Depolarizing(i; p = noise_p))
        push!(circuit, AmplitudeDamping(i; γ = noise_gamma))
        push!(circuit, PhaseDamping(i; γ = noise_gamma_phase))
    end
end

# -----------------------------
    # 3. Loss
    # -----------------------------
    loss_obj = LossExpectationRealDM(H, ψ0)
    θ = active_parameters(circuit)

    # -----------------------------
    # 4. Gradient-descent loop
    # -----------------------------
    η      = 0.01
    epochs = 1

    #logprintln(io, "Start training: L=$L, depth=$depth, η=$η, epochs=$epochs")
    #logprintln(io, "Log file: $(abspath(log_path))")

    for epoch in 1:epochs
        gθ = gradientMJVQC(loss_obj, circuit)
        gθ_ref = [ 7.81523780e-18, -4.22299617e-01,  3.29406982e-02,  7.81523780e-18,
 -1.61862269e-01,  1.27115982e-02, -1.00423504e-16,  2.57550759e-01,
  4.04886744e-02,  3.29406982e-02, -3.92636658e-01, -4.99577485e-03,
  2.16386042e-02, -3.49062349e-01, -5.45127770e-02,  4.50059994e-02,
 -1.73069409e-01,  8.98257086e-03, -4.99577485e-03, -3.24104369e-01,
  3.65028605e-03, -4.22764260e-02, -3.21447173e-01,  1.14307498e-02,
  5.43280702e-02, -1.70256168e-01, -1.50810359e-02]
        println("Computed gθ = ", gθ)
        println("Reference gθ = ", gθ_ref)

        @test length(gθ) == length(gθ_ref)
        @test isapprox(gθ, gθ_ref; rtol=1e-8, atol=1e-8)
        #logprintln(io, "Epoch $epoch: gradient = $gθ")
        #θ .-= η .* gθ
        #reset_parameters!(circuit, θ)

        #current_loss = loss_obj(circuit)
        #logprintln(io, "Epoch $epoch: loss = $current_loss")
        #logprintln(io, "")
    end

    #logprintln(io, "Finished training.")

end

