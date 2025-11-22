# docs/make.jl
using Documenter
using MyJuliVQC

makedocs(
    sitename = "MyJuliVQC.jl",
    modules  = [MyJuliVQC],
    format   = Documenter.HTML(;
        prettyurls = get(ENV, "CI", "false") == "true",
        canonical  = "https://HanDirac.github.io/MyJuliVQC.jl/",
    ),
    pages = [
        "Home" => "index.md",
        "Getting Started" => Any[
            "Installation" => "getting-started/installation.md",
            "Quick Start"  => "getting-started/quick-start.md",
        ],
        "Manual" => Any[
            "Initialize a quantum state"        => "manual/initialize-state.md",
            "Quantum gates"                     => "manual/quantum-gates.md",
            "Noise channels"                    => "manual/noise-channels.md",
            "Manipulating and running circuits" => "manual/circuits.md",
            "Qubit operators"                   => "manual/qubit-operators.md",
            "Automatic differentiation"         => "manual/automatic-differentiation.md",
            "Threading control (MyJuliVQC)"     => "manual/threading-control.md",
        ],
        "Examples" => Any[
            "Creating a Variational Quantum Circuit" => "examples/vqc-example.md",
            "A Simple Application with Flux"         => "examples/flux-example.md",
            "Utility Functions"                      => "examples/utilities.md",
        ],
    ],
    checkdocs = :none,   # ✅ 注意：在 pages 外面，这里是 makedocs 的关键字参数
)

deploydocs(
    repo      = "github.com/HanDirac/MyJuliVQC.jl.git",
    devbranch = "main",
)
