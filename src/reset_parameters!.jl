using LinearAlgebra
using ..MyJuliVQC: QuantumGate, QuantumMap, QCircuit

"""
        reset_parameters!(circ::QCircuit, paras::AbstractVector)

    Reset all trainable parameters of parametric operations in `circ` using
    the given scalar array `paras`, whose order must match that returned by
    `active_parameters(circ)`.

    After reset, those parametric operations are updated in-place (their
    internal `params` field is changed), and can later be instantiated by
    their `build` function into concrete quantum gates (or channels).

    - If the number of provided parameters is insufficient or excessive,
      an `ArgumentError` is thrown.
    - This function mutates `circ` in place (note the `!`) and returns
      the modified circuit itself.
    """
    function reset_parameters!(circ::QCircuit, paras::AbstractVector)
        # Use a mutable index ref to consume `paras` progressively during traversal
        idx = Ref(1)
        _reset_params_in_vec!(circ.operations, paras, idx)
        # If there are leftover parameters, throw an error instead of silently ignoring them
        if idx[] != length(paras) + 1
            needed = idx[] - 1
            throw(ArgumentError("too many parameters: expected $needed, got $(length(paras))"))
        end
        return circ
    end

    # Recursively process a vector of operations; if a nested QCircuit existed
    # we would descend into it; when encountering a parametric operation, we
    # update its parameters and keep its mask/build.
    function _reset_params_in_vec!(ops::Vector{CircuitElement}, paras::AbstractVector, idxref::Base.RefValue{Int})
        for i in eachindex(ops)
            op = ops[i]
            #if op isa QCircuit
                # Recurse into subcircuit
            #    _reset_params_in_vec!(op.operations, paras, idxref)
            if op isa ParamOp
                # Get original parameters and mask
                p0   = collect(op.params)          # copy of original params (preserving element type)
                mask = BitVector(op.mask)
                nvar = count(mask)

                # Check whether there are enough remaining parameters
                remain = length(paras) - idxref[] + 1
                if remain < nvar
                    throw(ArgumentError("not enough parameters: need $nvar more starting at position $(idxref[]), only $remain provided"))
                end

                # Build the parameter vector to be passed to `build`, overwriting
                # only the entries where mask is true
                p_use = similar(p0)
                copyto!(p_use, p0)
                idxs = findall(mask)
                @inbounds for k in eachindex(idxs)
                    j = idxs[k]
                    p_use[j] = convert(eltype(p0), paras[idxref[] + k - 1])
                end
                # Consume the number of parameters used in this step
                idxref[] += nvar

                # Update only `params`; keep mask and build unchanged
                ops[i] = ParamOp{nothing, eltype(p_use), typeof(op.build)}(p_use, mask, op.build)
            else
                # Fixed gates / quantum channels / other non-parametric objects: do nothing
            end
        end
        return nothing
    end