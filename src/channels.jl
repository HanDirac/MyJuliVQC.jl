module Channels

using LinearAlgebra
using ..MyJuliVQC: QuantumMap

export AmplitudeDamping, PhaseDamping, Depolarizing


function AmplitudeDamping(pos::Integer; γ,
        elty::Type=ComplexF64, check_tp::Bool=true)
    0 ≤ γ ≤ 1 || throw(ArgumentError("γ must be in [0,1], got $γ"))
    T = elty
    K0 = T[1 0; 0 sqrt(T(1)-T(γ))]
    K1 = T[0 sqrt(T(γ)); 0 0]
    QuantumMap(pos, (K0, K1); elty=T, check_tp=check_tp)
end

function PhaseDamping(pos::Integer; γ,
        elty::Type=ComplexF64, check_tp::Bool=true)
    0 ≤ γ ≤ 1 || throw(ArgumentError("γ must be in [0,1], got $γ"))
    T = elty
    K0 = T[1 0; 0 sqrt(T(1)-T(γ))]
    K1 = T[0 0; 0 sqrt(T(γ))]
    QuantumMap(pos, (K0, K1); elty=T, check_tp=check_tp)
end

function Depolarizing(pos::Integer; p,
        elty::Type=ComplexF64, check_tp::Bool=true)
    0 ≤ p ≤ 4/3 || throw(ArgumentError("p must be in [0,4/3], got $p"))
    T  = elty
    I2 = Matrix{T}(I, 2, 2)
    X  = T[0 1; 1 0]
    Y  = T[0 -im; im 0]
    Z  = T[1 0; 0 -1]
    a0 = sqrt(T(1) - T(3)*T(p)/T(4))
    a  = sqrt(T(p)/T(4))
    K0, K1, K2, K3 = a0 .* I2, a .* X, a .* Y, a .* Z
    QuantumMap(pos, (K0, K1, K2, K3); elty=T, check_tp=check_tp)
end

end # module
