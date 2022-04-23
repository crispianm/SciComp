module DiffEqToolkit

# Packages
using PlotlyJS
using NLsolve
using ForwardDiff
using Test
using LinearAlgebra
using ProgressMeter

# Subsystems
include("./examples/example_functions.jl")
include("ode_solver.jl")
include("numerical_shooting.jl")
include("numerical_continuation.jl")
include("finite_difference.jl")
include("visualisation.jl")

end