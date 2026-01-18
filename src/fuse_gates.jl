using LinearAlgebra
using ..MyJuliVQC: QuantumGate, QuantumMap, QCircuit

"""
        fuse_gates(circ::QCircuit) -> QCircuit

    Perform a **preliminary** simplification on the given circuit.  
    Only supports **noiseless** circuits (does not allow `QuantumMap`).  

    If there exists a single-qubit gate adjacent to a two-qubit gate
    and acting on one of the qubits of that two-qubit gate,
    the single-qubit gate will be **absorbed** into the two-qubit gate.  
    
    Only **adjacent** single-qubit gates are absorbed; the scan proceeds
    from left to right, and earlier two-qubit gates have higher absorption priority.  
    
    After simplification, all gates will be **non-parametric** `QuantumGate`s
    (if a gate was originally a specialized type, it will degrade into a general `QuantumGate` after fusion).
    """
    function fuse_gates(circ::QCircuit)
        # Work on a flattened copy to avoid in-place modification of the original circuit
        ops0 = circ

        # Noiseless check + materialize parametric operations into concrete QuantumGates
        materialized = Vector{Any}(undef, length(ops0))
        for (i, x) in enumerate(ops0)
            if x isa QuantumMap
                throw(ArgumentError("fuse_gates only supports noiseless circuits (found QuantumMap at position $i)"))
            elseif x isa ParamOp
                materialized[i] = x.build(x.params)  # Instantiate
            else
                materialized[i] = x
            end
            # Safety check: still disallow non-gate objects
            if !(materialized[i] isa QuantumGate)
                throw(ArgumentError("fuse_gates expects only gate-like ops after materialization; got $(typeof(materialized[i])) at #$i"))
            end
        end

        # Utility functions: determine 1-/2-qubit gates
        _arity(g::QuantumGate) = length(qubits(g))
        _is1(g::QuantumGate) = _arity(g) == 1
        _is2(g::QuantumGate) = _arity(g) == 2

        # Utility: create a 2×2 identity matrix of a given element type
        _I2(::Type{T}) where {T} = Matrix{T}(I, 2, 2)

        # Embed a single-qubit gate U into the two-qubit gate subspace (p1, p2):
        # side = :left  => single-qubit gate comes after the two-qubit gate: left multiply ((U⊗I) * T or (I⊗U) * T)
        # side = :right => single-qubit gate comes before the two-qubit gate: right multiply (T * (U⊗I) or T * (I⊗U))

        function _absorb!(M::Matrix{T}, p12::NTuple{2,Int}, U::Matrix{T}, q::Int; side::Symbol) where {T}
            p1, p2 = p12
            I2 = _I2(T)
            if q == p1
                K = kron(I2, U)  # U on faster slot
            elseif q == p2
                K = kron(U, I2)  # U on slower slot
            else
                return false  # Different qubit — cannot absorb
            end
            if side === :left
                Mul = K * M
            elseif side === :right
                Mul = M * K
            else
                throw(ArgumentError("unknown side=$side"))
            end
            # Overwrite with the new fused matrix
            M .= Mul
            return true
        end

        # For each two-qubit gate, first absorb the left neighbor (earlier in sequence)
        # if it acts on one of its qubits, then absorb the right neighbor (later in sequence)
        out = Vector{Any}()
        i = 1
        while i <= length(materialized)
            op = materialized[i]::QuantumGate
            if _is2(op)
                # Attempt to absorb single-qubit gates adjacent to this two-qubit gate
                M2  = Matrix{eltype(matrix(op))}(matrix(op))  # Copy to mutable matrix
                p12 = (qubits(op)[1], qubits(op)[2])

                # Absorb from the left (previous gate)
                if !isempty(out) && (out[end] isa QuantumGate) && _is1(out[end])
                    g1 = out[end]::QuantumGate
                    q1 = qubits(g1)[1]
                    # Absorb it (single-qubit gate comes “before”, so right multiply)
                    ok = _absorb!(M2, p12, matrix(g1), q1; side=:right)
                    if ok
                        pop!(out)  # Left neighbor absorbed
                    end
                end

                # Absorb from the right (next gate)
                if i < length(materialized)
                    nxt = materialized[i+1]
                    if (nxt isa QuantumGate) && _is1(nxt)
                        g1 = nxt::QuantumGate
                        q1 = qubits(g1)[1]
                        # Absorb it (single-qubit gate comes “after”, so left multiply)
                        ok = _absorb!(M2, p12, matrix(g1), q1; side=:left)
                        if ok
                            i += 1  # Skip this right neighbor
                        end
                    end
                end

                 # Output the fused two-qubit general gate (no longer a specialized type)
                push!(out, QuantumGate(qubits(op), M2; validate=false))
            else
                # Non-two-qubit gates (1-qubit or higher) — output directly
                push!(out, op)
            end
            i += 1
        end

        # Return a new QCircuit (containing only concrete QuantumGates)
        return QCircuit(out)
    end