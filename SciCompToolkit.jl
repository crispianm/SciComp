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
include("./src/FiniteDifference.jl")
include("./src/NumericalContinuation.jl")
include("./src/NumericalShooting.jl")
include("./src/ODESolver.jl")
include("./src/Visualisation.jl")

# Export all functions
@exportAll()

end
