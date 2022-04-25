module SciCompToolkit

# Packages
using ExportAll
using ForwardDiff
using LinearAlgebra
using NLsolve
using PlotlyJS
using ProgressMeter
using Test


# Subsystems
include("./examples/example_functions.jl")
include("finite_difference.jl")
include("numerical_continuation.jl")
include("numerical_shooting.jl")
include("ode_solver.jl")
include("visualisation.jl")

# Export all functions
@exportAll()

end
