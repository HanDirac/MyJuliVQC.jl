# Installation

## 1: Create a Dedicated Julia Environment (Highly Recommended)

We strongly recommend creating a dedicated Julia environment for MyJuliVQC.jl to prevent dependency conflicts with other projects. Follow these steps:

1. Open the Julia REPL.  
2. Navigate to your project folder (or any desired directory).  
3. Activate a new environment:  
   ```julia
   julia> ] activate .
   julia> ] instantiate
   ```
   This creates a new `Project.toml` and an empty environment in the current folder.

---

## 2. Install MyJuliVQC

MyJuliVQC can be installed directly from GitHub:
    ```julia
    julia> ] add https://github.com/HanDirac/MyJuliVQC.jl
    ```
Then load the package:
    ```julia
    julia> using MyJuliVQC
    ```
If this succeeds without errors, the installation is complete.

---





