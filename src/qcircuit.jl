using LinearAlgebra
using ..MyJuliVQC: QuantumGate, QuantumMap

# Determine whether something is a “parametric operation descriptor”
#_is_paramop(x) = x isa NamedTuple &&
#                all(k -> haskey(x, k), (:params, :mask, :build)) &&
#                (x.build isa Function)

# Determine whether an operation is allowed to be put into a circuit
_is_allowed_op(x) = (x isa QuantumGate) || (x isa QuantumMap) || (x isa QCircuit) || (x isa ParamOp)

# Flatten: turn nested QCircuit structures into a flat sequence of elements
#function _flatten_ops!(out::Vector{Any}, item)
#    if item isa QCircuit
#        for sub in item.ops
#            _flatten_ops!(out, sub)
#        end
#    else
#        push!(out, item)
#    end
#end

const CircuitElement = Union{QuantumGate, QuantumMap, ParamOp}

struct QCircuit <: QuantumOperation
    operations::Vector{CircuitElement}
end

QCircuit() = QCircuit(CircuitElement[])

# Accept “possibly nested” inputs; flatten upon construction
function QCircuit(xs::AbstractVector; flatten::Bool=true)
    if !flatten
        error("To keep types stable, it is recommended to always use flatten=true")
    end
    out = CircuitElement[]
    for x in xs
        _flatten_push!(out, x)
    end
    return QCircuit(out)
end

# Support pushing a mix of QCircuit/gates/ParamOp; flatten QCircuit immediately when encountered
function _flatten_push!(out::Vector{CircuitElement}, x)
    if x isa QCircuit
        for y in x.operations
            _flatten_push!(out, y)
        end
    elseif x isa CircuitElement
        push!(out, x)
    else
        throw(ArgumentError("Unsupported element: $(typeof(x))"))
    end
end

Base.push!(c::QCircuit, x) = (_flatten_push!(c.operations, x); c)

# Convenience constructors: empty circuit / single op / vararg form / iteration overload
#QCircuit() = QCircuit(Any[]; flatten=true, validate=false)
#QCircuit(op; kwargs...) = QCircuit([op]; kwargs...)
#QCircuit(op1, op2, ops...; kwargs...) = QCircuit(Any[op1, op2, ops...]; kwargs...)
Base.iterate(c::QCircuit, st...) = iterate(c.operations, st...)


# Common interfaces (container-like; do not execute anything)
Base.length(c::QCircuit) = length(c.operations)
Base.getindex(c::QCircuit, i::Int) = c.operations[i]

# Helper: obtain a **flattened** list of elements (does not modify the original object)
#flattened_ops(c::QCircuit) = (tmp = Any[]; _flatten_ops!(tmp, c); tmp)