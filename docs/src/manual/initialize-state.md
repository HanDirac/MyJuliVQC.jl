# Initialize State
The first step in using **MyJuliVQC** for quantum circuit simulation is to initialize a quantum state represented as a state vector. MyJuliVQC provides two constructors: `StateVector` for **pure states** and `DensityMatrix` for **mixed states**. 

From a mathematical viewpoint:
- An $n$-qubit **pure state** should be understood as a a vector of size $2^n$, corresponding to a rank-$n$ tensor with each index of dimension 2.
- An $n$-qubit **mixed state** is represented as a density matrix of size $2^n×2^n$, which can be viewed as a rank-$2n$ tensor where each dimension also has size 2.

Internal index conventions:
- For **pure states**, qubits are internally labeled from `1` to `n`.
- For **mixed states**, the _ket_ indices are labeled from `1` to `n`, and the corresponding _bra_ indices from `n+1` to `2n`.

Both state representations use **column-major storage**, meaning lower-index dimensions vary fastest.

These implementation details are generally irrelevant to typical users unless direct access to the raw data is required.

```julia
using MyJuliVQC

# Initialize |00⟩ for 2 qubits
state = StateVector(2)

n = 2
pure_state  = StateVector(n)
mixed_state = DensityMatrix(n)

# Custom pure state (user-specified vector)
custom_pure_state = StateVector([0.0, 0.1, 0.0, 0.0])

# Custom mixed state (user-specified flattened matrix)
# NOTE:
#   If a vector is used to construct a DensityMatrix,
#   **the vector must follow column-major (column-first) order**.
#   That is, the matrix is flattened by stacking all columns in sequence.
custom_mixed_state = DensityMatrix([0.5, -0.5im, 0.5im, 0.5]) #column-major order)

# -------------------------------------------------------
# Construct an equivalent mixed state using a 2×2 matrix
# -------------------------------------------------------

# 2×2 density matrix explicitly written
ρmat = [
    0.5      0.5im;
   -0.5im    0.5
]

# Construct DensityMatrix from matrix form
mixed_state_from_matrix = DensityMatrix(ρmat)

# -------------------------------------------------------
# Test equivalence
# -------------------------------------------------------

println("custom_mixed_state data = ", custom_mixed_state.data)
println("mixed_state_from_matrix data = ", mixed_state_from_matrix.data)

# They should match element-wise
println("Are the two density matrices equal? ",
        custom_mixed_state.data ≈ mixed_state_from_matrix.data)

```