# Code adapted from QuantumCircuits.jl by Weiyou Liao
# (the code written by the original author of the corresponding paper)
#
# The original paper does not clearly explain the implementation details
# of `QubitsTerm`, so we directly follow Weiyou Liao's implementation.
#
# QubitsTerm represents a single Pauli-string term:
#    positions::Vector{Int}          -- qubit indices (1-based), stored in ascending order
#    op::Vector{<:AbstractMatrix}    -- the 2×2 local operator on each qubit (aligned with positions)
#    coeff::Number                   -- scalar prefactor (coefficient), real or complex
#
# This file provides multiple convenient constructors (vector / tuple / dict / pairs)
# and common operations:
# - copy, isempty, adjoint (Hermitian conjugate), scalar multiplication, eltype
# - internal helpers: _get_normal_order, _get_op, dict_to_site_ops
#
# Design considerations:
# - (positions, ops) are normalized to ascending qubit order to ensure a canonical
#   representation (avoid equivalent representations with different ordering).
# - Strings such as "X", "Y", "Z", "+", "-" are supported and mapped to matrices
#   via the internal `_op_mapping` table.
# - `coeff` may be specified via a keyword argument, default = 1.0.
#
# Usage examples (in Julia REPL):
# ```
# julia> t1 = QubitsTerm([2,1], ["X","Z"], 0.5)   # positions and ops must have same length
# julia> t2 = QubitsTerm(1=>"X", 3=>"Y"; coeff=2.0)
# julia> t3 = QubitsTerm(Dict(2=>"Z", 1=>"X"))    # coeff defaults to 1.0
# julia> copy(t1)
# julia> adjoint(t2)
# julia> 3 * t3
# ```
#

struct QubitsTerm <: QuantumHamiltonian
	positions::Vector{Int}
	op::Vector{<:AbstractMatrix}
	coeff::Number

	function QubitsTerm(pos::Vector{Int}, m::Vector, v::Number)
		(length(pos) == length(m)) || error("number of sites mismatch with number of ops.")
		pos, m = _get_normal_order(pos, m)
		new(pos, m, v)
	end 
end

# -----------------------------
# Accessors (simple convenience helpers)
# -----------------------------
positions(x::QubitsTerm) = x.positions
oplist(x::QubitsTerm) = x.op
coeff(x::QubitsTerm) = x.coeff

# -----------------------------
# Convenience constructors
# -----------------------------
# Support tuple input for positions (internally converted to Vector)
QubitsTerm(pos::Tuple, m::Vector, v::Number) = QubitsTerm([pos...], m, v)

"""
    QubitsTerm(x::AbstractDict{Int}; coeff::Number=1.0)

Construct a QubitsTerm from a dictionary.
- Keys are qubit indices (`Int`).
- Values may be either:
    * a string operator ("X", "Z", ...)  
    * or a 2×2 matrix.
- `coeff` is an optional keyword argument (default 1.0).
Example:
    QubitsTerm(Dict(1=>"X", 3=>"Y"); coeff=2.0)
"""
function QubitsTerm(x::AbstractDict{Int}; coeff::Number=1.)
	sites, ops = dict_to_site_ops(x)
	return QubitsTerm(sites, ops, coeff)
end

"""
    QubitsTerm(i1 => op1, i2 => op2, ...; coeff=1.0)

Construct a `QubitsTerm` from several `(index => operator)` pairs.

Parameters
----------
- Each pair’s key must be an `Int` (1-based qubit index: `1`, `2`, ...).
- Each value may be:

  • **String** (`AbstractString`)  
    Must be a supported operator label, such as `"X"`, `"Y"`, `"Z"`, `"I"`,
    `"+"`, `"-"`, `"u"`, `"d"`, `"0"`, `"1"`.  
    These will be mapped to 2×2 matrices by `_op_mapping`.  
    **Case sensitive.**

  • **Matrix** (`AbstractMatrix`)  
    Any user-provided 2×2 operator. (Typical use: single-qubit operator matrices.)

- `...` indicates arbitrary many `(i => op)` pairs.

Keyword arguments
-----------------
- `coeff::Number = 1.0`  
  The numerical prefactor. Can be real or complex.

Notes
-----
- All pairs are collected into a `Dict` first; then construction proceeds
  through `QubitsTerm(Dict(...); coeff=...)`. This ensures normalization
  into ascending qubit order.
- If duplicate keys appear, the last one wins (Dict semantics).  
  Avoid specifying the same qubit more than once.
- If you intend to specify multi-qubit operators or matrices not of size 2×2,
  ensure your downstream usage supports it, as this structure is designed for
  Pauli-string–style one-qubit operators.

Examples
--------
```julia
# Using string operators
t1 = QubitsTerm(1=>"X", 3=>"Y"; coeff=0.5)

# Using explicit matrices
σx = [0. 1.; 1. 0.]
t2 = QubitsTerm(2=>σx, 4=>"Z"; coeff=1.0+0.0im)

# Equivalent: constructing from a Dict
t3 = QubitsTerm(Dict(1=>"X", 3=>"Y"); coeff=0.5)
"""
QubitsTerm(x::Pair{Int, <:Union{AbstractString, AbstractMatrix}}...; coeff::Number=1.) = QubitsTerm(
	Dict(x...), coeff=coeff)


Base.copy(x::QubitsTerm) = QubitsTerm(copy(positions(x)), copy(oplist(x)), coeff(x))
Base.isempty(x::QubitsTerm) = isempty(oplist(x))
Base.adjoint(x::QubitsTerm) = QubitsTerm(positions(x), [item' for item in oplist(x)], conj(coeff(x)))
Base.:*(x::QubitsTerm, y::Number) = QubitsTerm(positions(x), oplist(x), coeff(x)*y)
Base.:*(x::Number, y::QubitsTerm) = y * x

function Base.eltype(x::QubitsTerm)
	T = typeof(coeff(x))
	for item in oplist(x)
		T = promote_type(T, eltype(item))
	end
	return T
end

function _get_normal_order(key::Vector{Int}, op)
	seq = sortperm(key)
	return key[seq], [_get_op(item) for item in op[seq]]
end


_get_op(m::AbstractMatrix) = m
_get_op(m::String) = _op_mapping[m]

const _op_mapping = Dict("X"=>[0. 1.; 1. 0.], "Y" => [0. -im; im 0.], "Z"=>[1. 0.; 0. -1.], "+"=>[0. 1.; 0. 0.], "-"=>[0. 0.; 1. 0.], 
	"I"=>[1. 0.; 0. 1.], "u"=>[1. 0.; 0. 0.], "d"=>[0. 0.; 0. 1.], "0"=>[1. 0.; 0. 0.], "1"=>[0. 0.; 0. 1.])


function dict_to_site_ops(opstr::AbstractDict)
	sites = []
	ops = []
	for (k, v) in opstr
	    push!(sites, k)
	    push!(ops, v)
	end
	return [sites...], [ops...]
end