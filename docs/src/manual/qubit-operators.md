# Qubit Operators

In MyJuliVQC, multi-qubit observables and Hamiltonians are represented using two
closely related types:

- `QubitsTerm` – a single Pauli-string–style term (or more generally, a product of
  one-qubit operators with a scalar coefficient).
- `QubitsOperator` – a container holding a sum of `QubitsTerm` objects.

This design provides a flexible way to build
spin Hamiltonians, tensor-product observables, and general operators needed in
variational quantum algorithms.

Once a `QubitsOperator` `op` has been constructed, its expectation value on a state
(either `StateVector` or `DensityMatrix`) can be computed with:

```julia
expectation(op, state)
```

(implemented in a separate file; not in `qubits_term.jl` / `qubitsoperator.jl` directly).

---

## QubitsTerm: Single-Term Representation

A `QubitsTerm` represents a single product of local operators acting on specific qubits.

### Accessors

There are basic operations:

```julia
copy(t::QubitsTerm)
isempty(t::QubitsTerm)
adjoint(t::QubitsTerm)        # Hermitian conjugate (†), conjugating coeff and local matrices
3 * t                         # scalar multiplication (left or right)
Base.eltype(t::QubitsTerm)    # promoted numeric element type
```

---

## Constructing QubitsTerm

Several convenience constructors are provided.

### 1. From explicit vectors

```julia
QubitsTerm(pos::Vector{Int}, m::Vector, v::Number)
QubitsTerm(pos::Tuple,       m::Vector, v::Number)
```

Example:

```julia
using MyJuliVQC

pos = [2, 1]
ops = ["X", "Z"]    # will be mapped to matrices via _op_mapping
t   = QubitsTerm(pos, ops, 0.5)
```

Internally, `pos` and `ops` are normalized to ascending qubit order using an internal
helper `_get_normal_order`, so the canonical representation is on qubits `[1, 2]`.

### 2. From a dictionary

```julia
QubitsTerm(x::AbstractDict{Int}; coeff::Number = 1.0)
```

- keys: `Int` (1-based qubit indices),
- values: either strings (`"X"`, `"Z"`, …) or explicit `2×2` matrices.

Example:

```julia
using MyJuliVQC

t = QubitsTerm(Dict(1 => "X", 3 => "Y"); coeff = 2.0)
```

Internally, the dictionary is decomposed by `dict_to_site_ops` into
`(sites, ops)` and passed to the main constructor.

### 3. From pairs

This is often the most convenient syntax in user code:

```julia
QubitsTerm(i1 => op1, i2 => op2, ...; coeff::Number = 1.0)
```

Each pair uses:

- key: `Int` for the qubit index.
- value: either
  - an `AbstractString` (e.g. `"X"`, `"Y"`, `"Z"`, `"I"`, `"+"`, `"-"`,
    `"u"`, `"d"`, `"0"`, `"1"`), or
  - an `AbstractMatrix` (typically a `2×2` matrix).

These symbolic labels are mapped to matrices using the internal `_op_mapping` table:

```julia
const _op_mapping = Dict(
    "X" => [0. 1.; 1. 0.],
    "Y" => [0. -im; im 0.],
    "Z" => [1. 0.; 0. -1.],
    "+" => [0. 1.; 0. 0.],
    "-" => [0. 0.; 1. 0.],
    "I" => [1. 0.; 0. 1.],
    "u" => [1. 0.; 0. 0.],
    "d" => [0. 0.; 0. 1.],
    "0" => [1. 0.; 0. 0.],
    "1" => [0. 0.; 0. 1.],
)
```

Examples:

```julia
using MyJuliVQC

# Using string operators
t1 = QubitsTerm(1 => "X", 3 => "Y"; coeff = 0.5)

# Using explicit matrices
σx = [0. 1.; 1. 0.]
t2 = QubitsTerm(2 => σx, 4 => "Z"; coeff = 1.0 + 0.0im)

# Equivalent to constructing from Dict
t3 = QubitsTerm(Dict(1 => "X", 3 => "Y"); coeff = 0.5)
```

Notes:

- All pairs are collected into a `Dict` first and then normalized to ascending qubit order.
- If duplicate indices appear, the last one wins (dictionary semantics).

---

## QubitsOperator: Sum of Terms

`QubitsOperator` is a container for a **sum of QubitsTerm**, representing a Hamiltonian
or general observable.

### Constructors

```julia
QubitsOperator()                      # empty operator
QubitsOperator(t1, t2, ...)          # from varargs QubitsTerm
QubitsOperator(ts::Vector{<:QubitsTerm})
QubitsOperator(ts::AbstractVector)   # will be collected as Vector{QubitsTerm}
```

Example:

```julia
using MyJuliVQC

op = QubitsOperator()
push!(op, QubitsTerm(1 => "X", 2 => "Z"; coeff = 0.5))
```

### Collection-like interface

```julia
length(op::QubitsOperator)
isempty(op::QubitsOperator)
getindex(op::QubitsOperator, i::Int)
iterate(op::QubitsOperator)
copy(op::QubitsOperator)
push!(op::QubitsOperator, t::QubitsTerm)
```

You can also inspect it via `show`:

```julia
println(op)
# QubitsOperator with N term(s):
#  [1] positions = ..., coeff = ...
```

### Algebraic operations

```julia
adjoint(op::QubitsOperator)    # term-wise Hermitian conjugate
λ * op                         # scalar multiplication
op * λ                         # scalar multiplication (commutative)
A + B                          # concatenation of terms from A and B
```


---

## Simplification: `simplify!`

```julia
simplify!(op::QubitsOperator)
```

This in-place operation:

- merges terms having the **same** `positions` and the **same** sequence of local
  matrices (tested via an exact string-based key), and
- replaces them with a single term whose coefficient is the sum of the original
  coefficients.

Internally, it builds a `Dict{String,QubitsTerm}` keyed by a serialization of:

- the positions vector,
- and the flattened matrices in `oplist(t)`.

This is robust for canonical Pauli-string terms constructed via `"X"`, `"Y"`, etc.,
because their matrices are fixed. However:

> If the user supplies slightly different but mathematically equivalent matrices
> (e.g. `[0 1; 1 0]` vs `[0.0 1.0; 1.0 0.0]`), they may not be merged, since the
> string representations differ.

Typical usage:

```julia
using MyJuliVQC
op = QubitsOperator(
    QubitsTerm(1=>"X", 2=>"X"; coeff=1.0),
    QubitsTerm(1=>"X", 2=>"X"; coeff=2.0),
)

simplify!(op)
# Now op has a single term with coeff = 3.0
```

---

## Example: 1D Heisenberg Hamiltonian

The following function builds a standard 1D Heisenberg Hamiltonian on `L` qubits:

```julia
using MyJuliVQC

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
```

You can then use this operator with `expectation`:

```julia
ψ = StateVector(3)
H = heisenberg_1d(3)
E = expectation(H, ψ)
println("⟨ψ|H|ψ⟩ = ", E)
```

For mixed states:

```julia
ρ = DensityMatrix(3)
Eρ = expectation(H, ρ)
println("Tr(ρ H) = ", Eρ)
```

---

## Summary

- `QubitsTerm` encodes a single Pauli-string–like term with a scalar coefficient.
- `QubitsOperator` is a sum of such terms, suitable for Hamiltonians and observables.
- Constructors accept vectors, dictionaries, and `i => op` pairs, with string labels
  mapped to predefined 2×2 matrices.
- `simplify!` merges identical terms by summing coefficients.
- These structures integrate naturally with `expectation(op, state)` and the
  variational/gradient engine in MyJuliVQC.
