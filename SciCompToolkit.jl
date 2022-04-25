module SciCompToolkit

# Packages
using PlotlyJS
using NLsolve
using ForwardDiff
using Test
using LinearAlgebra
using ProgressMeter
using ExportAll
using Pkg

dependencies = [
    "ExportAll",
    "ForwardDiff",
    "IJulia",
    "LinearAlgebra",
    "NLsolve",
    "PlotlyJS",
    "ProgressMeter",
    "Test",
]

Pkg.add(dependencies)

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