# Threading Control (MyJuliVQC Extension)

MyJuliVQC provides additional controls for threading and numeric precision, beyond the original JuliVQC implementation.

```julia
using MyJuliVQC

MyJuliVQC.set_threading!(;
    outer_threads = true,
    dloc_threshold = 8,
    blas_threads = 1,
)
```

This page will describe:
- how outer threading is applied in gate application,
- interaction with BLAS threads,
- recommended settings for different hardware (laptop / workstation / HPC).
