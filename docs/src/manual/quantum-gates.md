# Quantum Gates

The second step of using **MyJuliVQC** is to build a quantum circuit, which requires specifying the elementary quantum gate operations (and, for noisy circuits, quantum channels).

The most general way to define a quantum gate in MyJuliVQC is via the constructor

`QuantumGate(positions, data)`

where:

- `positions` specifies the qubit indices that the gate acts on, e.g. `(1, 3)` for a two-qubit gate on qubits 1 and 3;
- `data` is the raw matrix of the gate, which should be a unitary matrix of size `2^k × 2^k` for a `k`-qubit gate.

Internally, MyJuliVQC uses **column-major** ordering (Julia’s default) for storing gate matrices and state data.


## Matrix Element Convention for `QuantumGate`

More specifically, suppose `QuantumGate(positions, data)` (where $\texttt{positions} = (j_1, j_2, \dots, j_l)$) represents an operator $\hat{O}$ acting on the Hilbert subspace associated with the qubits ${j_1}, {j_2}, \ldots, {j_l}$. Let
$|x_{j_1}=k_1,\space x_{j_2}=k_2,\space \cdots, \space x_{j_l}=k_l\rangle$
denote a computational basis vector in this Hilbert space. Then
$$\texttt{data}[a,b]
=
\langle
x_{j_1}=a_1,x_{j_2}=a_2,\dots,x_{j_l}=a_l
\mid
\hat{O}
\mid
x_{j_1}=b_1,x_{j_2}=b_2,\dots,x_{j_l}=b_l
\rangle,$$
where $a,b \in {1,2,\dots,2^l}$, and the binary expansions of $(a-1)$ and $(b-1)$ correspond to $(a_l,\dots,a_2,a_1)$ and $(b_l,\dots,b_2,b_1)$, respectively.

As a concrete example, consider a two-qubit controlled gate whose control qubit is the 2nd qubit and target qubit is the 4th qubit. If we take
`positions = (2,4)`
then the corresponding `data` matrix can be written as 
$$\texttt{data}
=
\begin{pmatrix}
1 & 0 & 0 & 0 \\
0 & u_{11} & 0 & u_{12} \\
0 & 0 & 1 & 0 \\
0 & u_{21} & 0 & u_{22}
\end{pmatrix},$$
where
$$\begin{pmatrix} u_{11} & u_{12} \\ u_{21} & u_{22} \end{pmatrix}$$
denotes the single-qubit unitary acting on the target qubit when the control qubit takes the value $1$.

If we swap the order of `positions`, i.e. take
`positions} = (4,2)`,
then `data` should instead be written as
$$\texttt{data}
=
\begin{pmatrix}
1 & 0 & 0 & 0 \\
0 & 1 & 0 & 0 \\
0 & 0 & u_{11} & u_{12} \\
0 & 0 & u_{21} & u_{22}
\end{pmatrix},$$

If we change the control qubit of this two-qubit controlled gate to be the 4th qubit and the target qubit to be the 2nd qubit, then we may set
$$\texttt{positions}=(4,2),\space
\texttt{data}
=
\begin{pmatrix}
1 & 0 & 0 & 0 \\
0 & u_{11} & 0 & u_{12} \\
0 & 0 & 1 & 0 \\
0 & u_{21} & 0 & u_{22}
\end{pmatrix},$$
or
$$\texttt{positions}=(2,4),\space
\texttt{data}
=
\begin{pmatrix}
1 & 0 & 0 & 0 \\
0 & 1 & 0 & 0 \\
0 & 0 & u_{11} & u_{12} \\
0 & 0 & u_{21} & u_{22}
\end{pmatrix}.$$

---

## Built-in Quantum Gates

For convenience, MyJuliVQC also provides named constructors for commonly used quantum gates:

`XGate`, `YGate`, `ZGate`, `SGate`, `HGate`, `sqrtXGate`, `sqrtYGate`, `TGate`, `SWAPGate`, `iSWAPGate`, `CZGate`, `CNOTGate`, `TOFFOLIGate`, `FREDKINGate`, `RxGate`, `RyGate`,`RzGate`, `CRxGate`, `CRyGate`, `CRzGate`, `FSIMGate`. 

In the current implementation of **MyJuliVQC**, these special gates are implemented as
thin wrappers around the generic `QuantumGate` constructor. This design keeps the code
transparent and easy to read. Future versions may introduce specialized, highly optimized
kernels for some of these gates, but the public API will remain the same.

---

## Controlled Gates

MyJuliVQC supports generic controlled operations via dedicated constructors:

- `CONTROLGate(i, j, data)`: a controlled single-qubit gate  
  - `i`: control qubit  
  - `j`: target qubit  
  - `data`: the `2×2` matrix of the target single-qubit operation

- `CONTROLCONTROLGate(i, j, k, data)`: a doubly controlled single-qubit gate  
  - `i`, `j`: control qubits  
  - `k`: target qubit  
  - `data`: again the `2×2` matrix acting on the target qubit

These constructors are convenient when you want to promote an arbitrary single-qubit gate to a (multi-)controlled version without manually building the full `2^n × 2^n` matrix.

---

## Parametric Quantum Gates

Many variational algorithms require **parameterized gates**.  
In MyJuliVQC, such gates follow a unified interface

`G(i..., paras; isparas)`

where:

- `G` is a gate constructor such as `RxGate`, `RyGate`, `RzGate`, `FSIMGate`, etc.;
- `i...` are the qubit indices;
- `paras` is either:
  - a single scalar (for one-parameter gates), or
  - an array of scalars (for multi-parameter gates);
- `isparas` is a Boolean keyword:
  - `isparas = false`: the gate is treated as a **fixed** numerical gate;
  - `isparas = true`: the gate is treated as having **optimizable parameters**, and its
    parameters will be tracked by the variational/gradient engine.

---

## Examples: Non-Parametric and Parametric Gates

The following code illustrates how to construct basic non-parametric and parametric gates
in MyJuliVQC:

```julia
using MyJuliVQC

# Single-qubit X gate on qubit 1
n = 1
X = XGate(n)

# Two-qubit CNOT: control = 1, target = 2
ncontrol = 1
ntarget  = 2
CNOT = CNOTGate(ncontrol, ntarget)

# Rx gate: non-parametric vs parametric
θ = π / 2

# a non-parametric Rx gate (angle fixed, not tracked as a variational parameter)
non_para_Rx = RxGate(n, θ; isparas = false)

# a parametric Rx gate (its angle is treated as a variational parameter)
para_Rx = RxGate(n, θ; isparas = true)
```

