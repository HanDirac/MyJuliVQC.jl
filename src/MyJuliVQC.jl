module MyJuliVQC

<<<<<<< HEAD
using LinearAlgebra
using Random
=======
using LinearAlgebra, TensorOperations
>>>>>>> 2a33fd454380f1bd988e284a5fedd1e28e2a6bb6

# Unified exports
export StateVector, DensityMatrix, QCircuit
export apply!, apply, measure!, expectation
export QubitsOperator, QuantumOperation, QuantumPrimitive, QuantumHamiltonian
export active_parameters, reset_parameters!, fuse_gates
export QuantumGate, QuantumMap, n_qubits_mat_from_external

# === Global threading configuration (for internal read / user settings) ===
const _USE_OUTER_THREADS        = Ref(true)  # Whether to enable outer-level parallelism (for gather/scatter using @threads)
const _DLOC_THREAD_THRESHOLD    = Ref(8)     # Small-gate threshold: outer threading is considered only if dloc <= this value

"User API: configure threading strategy of MyJuliVQC; optionally also set the BLAS thread count."
function set_threading!(; outer_threads::Bool=_USE_OUTER_THREADS[],
                           dloc_threshold::Int=_DLOC_THREAD_THRESHOLD[],
                           blas_threads::Union{Nothing,Int}=nothing)
    _USE_OUTER_THREADS[]     = outer_threads
    _DLOC_THREAD_THRESHOLD[] = dloc_threshold
    if blas_threads !== nothing
        try
            LinearAlgebra.BLAS.set_num_threads(blas_threads)
        catch
            # Some platforms lack BLAS thread control; ignore silently
        end
    end
    return nothing
end

# Accessors for internal use (not exported)
use_outer_threads()      = _USE_OUTER_THREADS[]
dloc_thread_threshold()  = _DLOC_THREAD_THRESHOLD[]

export set_threading!

####  StateVector ####
include("statevector.jl")
    
####  DensityMatrix ####
include("densitymatrix.jl")

#### QuantumOperation ####
include("quantum_operation.jl")

####  QuantumGate ####
include("quantumgate.jl")
    
####  n_qubits_mat_from_external #### (development postponed)

####  Load gate library ####
include("gates.jl")
using .Gates: ParamOp
using .Gates: XGate, YGate, ZGate, HGate, SGate, TGate, sqrtXGate, sqrtYGate,
        RxGate, RyGate, RzGate, CRxGate, CRyGate, CRzGate,
        SWAPGate, iSWAPGate, CNOTGate, CZGate,
        TOFFOLIGate, FREDKINGate,
        CONTROLGate, CONTROLCONTROLGate,
        FSIMGate

export XGate, YGate, ZGate, HGate, SGate, TGate, sqrtXGate, sqrtYGate,
    RxGate, RyGate, RzGate, CRxGate, CRyGate, CRzGate,
    SWAPGate, iSWAPGate, CNOTGate, CZGate,
    TOFFOLIGate, FREDKINGate,
    CONTROLGate, CONTROLCONTROLGate,
    FSIMGate

####  QuantumMap ####
include("quantummap.jl")

####  Load quantum channels ####
include("channels.jl")
using .Channels
export AmplitudeDamping, PhaseDamping, Depolarizing

####  QCircuit ####
include("qcircuit.jl")
export CircuitElement

####  active_parameters ####
include("active_parameters.jl")
    
####  reset_parameters! ####
include("reset_parameters!.jl")

####  fuse_gates ####
include("fuse_gates.jl")

####  apply ####
include("apply.jl")

####  measure! ####
include("measure!.jl")

####  QubitsTerm ####
include("qubits_term.jl")
export QubitsTerm

####  QubitsOperator ####
include("qubitsoperator.jl")

####  expectation ####
include("expectation.jl")

####  gradient ####
include("gradient.jl")
<<<<<<< HEAD
export gradient, LossExpectationRealSV, LossExpectationRealDM
=======
export gradient
>>>>>>> 2a33fd454380f1bd988e284a5fedd1e28e2a6bb6




    

end
