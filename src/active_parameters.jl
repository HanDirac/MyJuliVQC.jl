using LinearAlgebra
using ..MyJuliVQC: QuantumGate, QuantumMap, QCircuit
"""
        active_parameters(c::QCircuit) -> Vector{Float64}

    Extracts all **active** variational parameters from the circuit `c` in order.
    - Only collects "parameterized operations" of the form `(; params::Vector, mask::BitVector, build::Function)`;
    - A parameter is considered active only when the corresponding `mask[i] == true`;
    - The return order matches the order of operations in the circuit
      (if there are nested structures, parameters are flattened according to their traversal order).
"""
    function active_parameters(c::QCircuit)
        vals = Float64[]
        for op in c  
            if op isa ParamOp
                @inbounds for i in eachindex(op.params, op.mask)
                    if op.mask[i]
                        p = op.params[i]
                        # Normalize to Float64; if it's complex with ≈0 imaginary part, take the real part
                        if p isa Real
                            push!(vals, float(p))
                        else
                            push!(vals, float(real(p)))
                        end
                    end
                end
            end
        end
        return vals
    end