using LinearAlgebra

# src/qubitsoperator.jl
# QubitsOperator: a sum of QubitsTerm (Pauli-string style Hamiltonian).
# NOTE: this file deliberately does NOT implement `expectation` here.
#       Implement `expectation(op, state)` later in a separate file when ready.

"""
    QubitsOperator(terms...)

Container storing a summation of `QubitsTerm` objects.

Constructors:
- `QubitsOperator()` -> empty operator
- `QubitsOperator(t1, t2, ...)` -> construct from varargs QubitsTerm
- `QubitsOperator(ts::Vector{<:QubitsTerm})` -> construct from a vector of terms

Example:
```julia
op = QubitsOperator()
push!(op, QubitsTerm(1=>"X", 2=>"Z"; coeff=0.5))
"""
struct QubitsOperator <: QuantumHamiltonian
    terms::Vector{QubitsTerm}
    function QubitsOperator(ts::AbstractVector{<:QubitsTerm})
        new(QubitsTerm[ts...]) # force concrete Vector{QubitsTerm}
    end
end

# convenience constructors
QubitsOperator() = QubitsOperator(QubitsTerm[])
QubitsOperator(ts::QubitsTerm...) = QubitsOperator(collect(ts))
# Allow passing Vector{Any} / Vector{T} as long as the elements are QubitsTerm
QubitsOperator(ts::AbstractVector) = QubitsOperator(collect(QubitsTerm, ts))

# Basic collection-like interface
Base.length(op::QubitsOperator) = length(op.terms)
Base.isempty(op::QubitsOperator) = isempty(op.terms)
Base.getindex(op::QubitsOperator, i::Int) = op.terms[i]
Base.iterate(op::QubitsOperator, s=1) = s > length(op) ? nothing : (op.terms[s], s+1)
Base.copy(op::QubitsOperator) = QubitsOperator(copy(op.terms))
Base.push!(op::QubitsOperator, t::QubitsTerm) = (push!(op.terms, t); op)

# Algebraic convenience operations
Base.adjoint(op::QubitsOperator) = QubitsOperator([adjoint(t) for t in op.terms])
Base.:*(λ::Number, op::QubitsOperator) = QubitsOperator([QubitsTerm(positions(t), oplist(t), coeff(t)*λ) for t in op.terms])
Base.:*(op::QubitsOperator, λ::Number) = λ * op
Base.:+(A::QubitsOperator, B::QubitsOperator) = QubitsOperator(vcat(A.terms, B.terms))

function Base.show(io::IO, op::QubitsOperator)
println(io, "QubitsOperator with ", length(op), " term(s):")
for (k,t) in enumerate(op.terms)
println(io, " [$k] positions=", positions(t), ", coeff=", coeff(t))
end
end

# Iteration metadata: QubitsOperator iterates over QubitsTerm
Base.IteratorEltype(::Type{QubitsOperator}) = Base.HasEltype()
Base.IteratorSize(::Type{QubitsOperator})   = Base.HasLength()
Base.eltype(::Type{QubitsOperator}) = QubitsTerm
Base.eltype(::QubitsOperator)       = QubitsTerm

"""
    numeric_eltype(op::QubitsOperator)

Return the promoted scalar element type used in numerical computations
(combining coefficients and the eltypes of all local matrices).
"""
function numeric_eltype(op::QubitsOperator)
    if isempty(op)
        return ComplexF64
    end
    T = typeof(coeff(op.terms[1]))
    for t in op.terms
        T = promote_type(T, Base.eltype(t))  # Base.eltype(t) comes from QubitsTerm's
    end
    return T
end

# -----------------------------------------------------------------------------
# Simplify: merge terms that have identical positions & oplist (exact equality)
# -----------------------------------------------------------------------------
"""
simplify!(op::QubitsOperator)

In-place combine identical terms (same positions and same sequence of small matrices)
by summing their coefficients. Returns the (modified) operator.
Note: this performs an exact equality check on matrix elements (suitable for string-built Pauli terms)
（This is fine for Pauli-string terms created from predefined symbols
("X", "Y", "Z", etc.), since their matrices are canonical.
However, if the user manually supplies matrices (e.g. [0 1; 1 0]
versus [0.0 1.0; 1.0 0.0]), then mathematically they represent the same operator,
but simplify! may treat them as different because their textual
representations differ.）.
"""
function simplify!(op::QubitsOperator)
    if isempty(op)
        return op
    end
    d = Dict{String, QubitsTerm}()
    for t in op.terms
        mats = oplist(t)
        matvecs = join([string(vec(m)) for m in mats], ";")
        key = string(positions(t)) * "|" * matvecs
        if haskey(d, key)
            existing = d[key]
            d[key] = QubitsTerm(positions(existing), oplist(existing),
                                coeff(existing) + coeff(t))
        else
            d[key] = t
        end
    end
    # modify in-place rather than reassign the field
    empty!(op.terms)
    append!(op.terms, values(d))
    return op
end