using MyJuliVQC

const gradientMJVQC = MyJuliVQC.gradient

log_path = "grad_desc_log_noisy.txt"
io = open(log_path, "w")   

function logprintln(io::IO, args...)
    println(args...)           # stdout
    println(io, args...)       # file
    flush(io)                  
end

try

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
    epochs = 1000

    logprintln(io, "Start training: L=$L, depth=$depth, η=$η, epochs=$epochs")
    logprintln(io, "Log file: $(abspath(log_path))")

    for epoch in 1:epochs
        gθ = gradientMJVQC(loss_obj, circuit)
        logprintln(io, "Epoch $epoch: gradient = $gθ")
        θ .-= η .* gθ
        reset_parameters!(circuit, θ)

        current_loss = loss_obj(circuit)
        logprintln(io, "Epoch $epoch: loss = $current_loss")
        logprintln(io, "")
    end

    logprintln(io, "Finished training.")

finally
    close(io)  
end