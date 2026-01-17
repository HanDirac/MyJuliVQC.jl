# Noise Channels

In addition to unitary quantum gates, realistic quantum circuits must also account for **noise**.

Noise processes are mathematically described by **quantum channels**, which generalize unitary evolution and allow non-unitary dynamics such as decoherence and relaxation.

In MyJuliVQC, quantum channels are represented using the **Kraus operator formalism**, providing a flexible and unified way to model general noisy processes.

---

## Defining General Quantum Channels

MyJuliVQC offers a universal constructor:

```julia
QuantumMap(positions, kraus)
```

where:

- `positions` specifies the qubits on which the channel acts (e.g., `(1,)`, `(2, 3)`);

- `kraus` is a vector of matrices

$$
K_{\ell}$$

representing the channel through

$$
\rho \mapsto \sum_{\ell} K_{\ell}\,\rho\,K_{\ell}^{\dagger}.$$


This interface allows the user to define arbitrary custom noise models, from single-qubit decoherence to multi-qubit correlated channels.

Example (defining a custom single-qubit channel):

```julia
using MyJuliVQC

K1 = [1 0; 0 sqrt(0.8)]
K2 = [0 sqrt(0.2); 0 0]

my_channel = QuantumMap(1, [K1, K2])   
```
---

## Built-in Noise Channels

For convenience, MyJuliVQC provides several commonly used **single-qubit** quantum channels as simple wrappers of QuantumMap.
These include:

- `AmplitudeDamping(pos; γ)`
- `PhaseDamping(pos; γ)`
- `Depolarizing(pos; p)`

Each constructor automatically generates the corresponding Kraus operators and returns a `QuantumMap`.
(Internally, these functions do not use any hardware-specific optimizations; they are implemented in a clean and transparent way suitable for research and educational use.)

Example:
```julia
using MyJuliVQC

channel = AmplitudeDamping(1; γ=0.1)
channelb = PhaseDamping(2; γ=0.2)
channelc = Depolarizing(3; p=0.3)
```
