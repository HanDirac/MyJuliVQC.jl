# Threading Control (MyJuliVQC Extension)

MyJuliVQC provides an optional threading-control interface that gives users more flexibility than the original JuliVQC implementation. It allows fine-grained control over how parallelism is applied during quantum-state evolution, which is important for large circuits or multi-core/HPC environments.

The configuration is done through:

```
MyJuliVQC.set_threading!(;
    outer_threads = true,
    dloc_threshold = 8,
    blas_threads = nothing,
)
```

This document explains:
- what each threading parameter means,
- how MyJuliVQC applies threading internally,
- interaction with BLAS-level threading,
- recommended settings for laptops, workstations, and HPC clusters.

---

## Global Threading Strategy

MyJuliVQC uses two potential sources of parallelism:

1. **Outer-level threading**  
   When applying small gates (1–2 qubits) to a large state vector, MyJuliVQC may distribute the work across Julia threads. This is also used for some density-matrix kernels.

2. **BLAS-level threading**  
   Some operations, especially involving large density matrices, rely on BLAS routines (MKL, OpenBLAS, etc.), which may themselves use multi-threading internally.

Using both at once may cause oversubscription. The function `set_threading!` helps control this balance.

---

## Configuration Interface

### Parameter: outer_threads

Controls whether MyJuliVQC uses Julia-level `Threads.@threads` for gate-application loops.

```
outer_threads = true    # enable outer-level parallel loops
outer_threads = false   # disable outer-level parallelism
```

Situations where setting this to `false` is beneficial:
- when BLAS is using many threads,
- when oversubscription occurs,
- when benchmarking indicates BLAS dominates performance.

---

### Parameter: dloc_threshold

Controls when outer threading is allowed.

- `dloc` = local dimension of a gate = `2^k` for a `k`-qubit gate.
- If `dloc <= dloc_threshold`, MyJuliVQC may activate outer threading.

Example meanings:

```
dloc_threshold = 8    # allow threading for gates up to 3 qubits
dloc_threshold = 4    # allow threading only for 1–2 qubit gates
```

---

### Parameter: blas_threads

Optional argument to configure BLAS:

```
blas_threads = N       # attempt to set BLAS to use N threads
blas_threads = nothing # leave BLAS thread count unchanged (default)
```

Some BLAS installations may not support runtime thread changes; MyJuliVQC silently ignores errors in such cases.

Setting `blas_threads = 1` is often useful if relying mostly on Julia threads.

---

## Example Usage

```
using MyJuliVQC

MyJuliVQC.set_threading!(;
    outer_threads = true,
    dloc_threshold = 8,
    blas_threads = 1,
)
```

This configuration:
- enables outer-level threading,
- applies it to small gates,
- restricts BLAS to a single thread.

---

## Internal Usage

MyJuliVQC stores configuration using:

```
_USE_OUTER_THREADS[]        # reflective boolean flag
_DLOC_THREAD_THRESHOLD[]    # integer threshold for small gates
```

and provides internal helper functions:

```
use_outer_threads()
dloc_thread_threshold()
```

These are not exported, but used by the kernel implementation.

---

## Recommended Configurations

### Laptops (4–8 cores)

```
set_threading!(;
    outer_threads = true,
    dloc_threshold = 8,
    blas_threads = 1,
)
```

### Workstations (16–32 cores)

```
set_threading!(;
    outer_threads = true,
    dloc_threshold = 12,
    blas_threads = 4,
)
```

### HPC Nodes (32–128+ cores)

```
set_threading!(;
    outer_threads = false,
    dloc_threshold = 16,
    blas_threads = 16,
)
```

HPC nodes typically benefit from strong BLAS kernels rather than Julia-threaded outer loops.

---

## Notes

- MyJuliVQC does not yet include specialized gate kernels like the original JuliVQC, so threading strategy can significantly affect performance.
- Users should benchmark real workloads to determine optimal settings.
- Future versions may include adaptive heuristics for automatically balancing Julia-level and BLAS-level threading.

