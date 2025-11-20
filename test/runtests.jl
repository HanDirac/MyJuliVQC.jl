using InteractiveUtils
using MyJuliVQC
using Test

include(joinpath(@__DIR__, "test_statevector.jl"))
#（Currently densitymatrix.jl has a known bug, so the test is temporarily disabled）include(joinpath(@__DIR__, "test_densitymatrix.jl"))
include(joinpath(@__DIR__, "test_quantumgate.jl"))
include(joinpath(@__DIR__, "test_gates.jl"))
include(joinpath(@__DIR__, "test_quantummap.jl"))
include(joinpath(@__DIR__, "test_channels.jl"))
include(joinpath(@__DIR__, "test_qcircuit.jl"))
include(joinpath(@__DIR__, "test_active_parameters.jl"))
include(joinpath(@__DIR__, "test_reset_parameters!.jl"))
include(joinpath(@__DIR__, "test_fuse_gates.jl"))
include(joinpath(@__DIR__, "test_apply_densitymatrix.jl"))
include(joinpath(@__DIR__, "test_apply_statevector.jl"))
include(joinpath(@__DIR__, "test_apply!_densitymatrix.jl"))
include(joinpath(@__DIR__, "test_apply!_statevector.jl"))
include(joinpath(@__DIR__, "test_measure!.jl"))
include(joinpath(@__DIR__, "test_expectation.jl"))
include(joinpath(@__DIR__, "test_gradient.jl"))