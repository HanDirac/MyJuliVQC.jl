using LinearAlgebra
using ..MyJuliVQC: QuantumGate, QuantumMap, QCircuit, StateVector, DensityMatrix

"""
        measure!(ψ::StateVector, i::Integer; rng=Random.GLOBAL_RNG) -> (outcome::Int, prob::Float64)

    Perform a computational-basis measurement on the `i`-th qubit, sample an outcome
    according to the probabilities, and **collapse the pure state in place**.
    Returns (measurement outcome 0/1, exact probability of that outcome).

    Convention: qubit 1 is the global LSB (least significant bit).
"""
    
    function measure!(ψ::StateVector{T}, i::Integer; rng=nothing) where {T<:Number}
        rng = rng === nothing ? Random.GLOBAL_RNG : rng
        data = ψ.data
        N = length(data)
        @assert ispow2(N) "state length must be 2^n"
        n = floor(Int, log2(N))
        @assert 1 ≤ i ≤ n "qubit index $i out of range 1..$n"

        stride = 1 << (i - 1)         # stride associated with this qubit
        block  = stride << 1          # size of two adjacent blocks
        p1 = zero(real(T))            # total probability of outcome == 1

        # First compute p1
        @inbounds for base = 0:block:(N-1)
            # For this block:
            #   front half (bit_i=0): indices base + 0 .. base + stride - 1
            #   back half  (bit_i=1): indices base + stride .. base + block - 1
            for k = 0:(stride-1)
                amp = data[base + stride + k + 1]
                p1 += real(amp*conj(amp))
            end
        end
        p1 = Float64(p1)
        p0 = max(0.0, 1.0 - p1)       # numerically robust

        # Sample the outcome
        b = (rand(rng) < p1) ? 1 : 0
        p = (b == 1) ? p1 : p0
        if p == 0.0
            # Zero probability for the chosen branch implies the other branch must occur
            b = 1 - b
            if b == 1
                @inbounds for base = 0:block:(N-1)
                    # Zero out the bit_i=0 half block
                    for k = 0:(stride-1)
                        data[base + k + 1] = zero(T)
                    end
                    # Keep the bit_i=1 half block (no need to renormalize; usually already unit norm)
                end
            else
                @inbounds for base = 0:block:(N-1)
                    # Keep the bit_i=0 half block
                    # Zero out the bit_i=1 half block
                    for k = 0:(stride-1)
                        data[base + stride + k + 1] = zero(T)
                    end
                end
            end
            return (b, 1.0)  # Early return; standard normalization logic below is skipped
        end

        # Collapse and renormalize in place
        normfac = inv(sqrt(p))
        if b == 1
            @inbounds for base = 0:block:(N-1)
                # Zero out the bit_i=0 half block
                for k = 0:(stride-1)
                    data[base + k + 1] = zero(T)
                end
                # Renormalize the bit_i=1 half block
                for k = 0:(stride-1)
                    data[base + stride + k + 1] *= normfac
                end
            end
        else
            @inbounds for base = 0:block:(N-1)
                # Renormalize the bit_i=0 half block
                for k = 0:(stride-1)
                    data[base + k + 1] *= normfac
                end
                # Zero out the bit_i=1 half block
                for k = 0:(stride-1)
                    data[base + stride + k + 1] = zero(T)
                end
            end
        end
        return (b, p)
    end

    """
        measure!(ρ::DensityMatrix, i::Integer; rng=Random.GLOBAL_RNG) -> (outcome::Int, prob::Float64)

    Perform a computational-basis measurement on the `i`-th qubit, sample an outcome
    according to the probabilities, and **collapse the density matrix in place**.
    Returns (measurement outcome 0/1, exact probability of that outcome).

    The implementation treats `ρ` as a `2^n × 2^n` matrix and keeps only the diagonal
    block consistent with the measurement result:  
    `ρ ← Π_b ρ Π_b / p_b`，where `Π_b = |b⟩⟨b|` acts on the `i`-th qubit
    (and identity acts on all other qubits).
    """
    function measure!(ρ::DensityMatrix{T}, i::Integer; rng=nothing) where {T<:Number}
        rng = rng === nothing ? Random.GLOBAL_RNG : rng
        d = 1 << nqubits(ρ)                 # 2^n
        M = reshape(ρ.data, d, d)           # view sharing underlying storage
        @assert 1 ≤ i ≤ floor(Int, log2(d)) "qubit index out of range"

        stride = 1 << (i - 1)
        block  = stride << 1

        # Probability p1 = tr(Π1 * ρ) = ∑_{row bit_i=1} ρ[row,row]
        p1 = 0.0
        @inbounds for base = 0:block:(d-1)
            for k = 0:(stride-1)
                idx = base + stride + k + 1
                p1 += real(M[idx, idx])
            end
        end
        p1 = Float64(p1)

        p0 = max(0.0, 1.0 - p1)
        b  = (rand(rng) < p1) ? 1 : 0
        p  = (b == 1) ? p1 : p0
        if p == 0.0
            # Zero probability for the chosen branch implies the other branch must occur
            b = 1 - b
            #invp = 1.0  # normalization factor (not needed)
                if b == 1
                @inbounds for rbase = 0:block:(d-1)
                    for cbase = 0:block:(d-1)
                        # Four combinations: 00, 01, 10, 11 sub-blocks; keep only 11
                        # 00
                        for k = 0:(stride-1), l = 0:(stride-1)
                            M[rbase+k+1, cbase+l+1] = zero(T)
                        end
                        # 01
                        for k = 0:(stride-1), l = 0:(stride-1)
                            M[rbase+k+1, cbase+stride+l+1] = zero(T)
                        end
                        # 10
                        for k = 0:(stride-1), l = 0:(stride-1)
                            M[rbase+stride+k+1, cbase+l+1] = zero(T)
                        end
                        # 11 -> normalization
                        #for k = 0:(stride-1), l = 0:(stride-1)
                            #M[rbase+stride+k+1, cbase+stride+l+1] *= invp
                        #end
                    end
                end
            else
                @inbounds for rbase = 0:block:(d-1)
                    for cbase = 0:block:(d-1)
                        # 00 -> normalization
                        #for k = 0:(stride-1), l = 0:(stride-1)
                        #    M[rbase+k+1, cbase+l+1] *= invp
                        #end
                        
                        # 01 / 10 / 11 -> zero out
                        for k = 0:(stride-1), l = 0:(stride-1)
                            M[rbase+k+1, cbase+stride+l+1] = zero(T)
                            M[rbase+stride+k+1, cbase+l+1] = zero(T)
                            M[rbase+stride+k+1, cbase+stride+l+1] = zero(T)
                        end
                    end
                end
            end
            return (b, 1.0)
        end

        # ρ ← Π_b ρ Π_b / p ：keep only the submatrix with both row and column
        invp = 1.0 / p
        if b == 1
            @inbounds for rbase = 0:block:(d-1)
                for cbase = 0:block:(d-1)
                    # Four combinations: 00, 01, 10, 11 sub-blocks; keep only 11
                    # 00
                    for k = 0:(stride-1), l = 0:(stride-1)
                        M[rbase+k+1, cbase+l+1] = zero(T)
                    end
                    # 01
                    for k = 0:(stride-1), l = 0:(stride-1)
                        M[rbase+k+1, cbase+stride+l+1] = zero(T)
                    end
                    # 10
                    for k = 0:(stride-1), l = 0:(stride-1)
                        M[rbase+stride+k+1, cbase+l+1] = zero(T)
                    end
                    # 11 -> renormalize
                    for k = 0:(stride-1), l = 0:(stride-1)
                        M[rbase+stride+k+1, cbase+stride+l+1] *= invp
                    end
                end
            end
        else
            @inbounds for rbase = 0:block:(d-1)
                for cbase = 0:block:(d-1)
                    # 00 -> renormalize
                    for k = 0:(stride-1), l = 0:(stride-1)
                        M[rbase+k+1, cbase+l+1] *= invp
                    end
                    # 01 / 10 / 11 -> zero out
                    for k = 0:(stride-1), l = 0:(stride-1)
                        M[rbase+k+1, cbase+stride+l+1] = zero(T)
                        M[rbase+stride+k+1, cbase+l+1] = zero(T)
                        M[rbase+stride+k+1, cbase+stride+l+1] = zero(T)
                    end
                end
            end
        end
        return (b, p)
    end