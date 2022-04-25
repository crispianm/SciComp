# About
Software suite of various numerical ODE and PDE solution methods, written in Julia for the EMAT30008 (Scientific Computing) module.


## Installation
Using the latest version of Julia prompt, run "include("requirements.jl")", which will install the following dependencies:
* ExportAll
* ForwardDiff
* LinearAlgebra
* NLsolve
* PlotlyJS
* ProgressMeter
* Test

Then, run SciCompToolkit.jl to install the module.

## Usage
After installation, use the module by calling "using .SciCompToolkit" in your Julia file.
Note: This does not work in IJulia notebooks for some reason.