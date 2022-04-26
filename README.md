## About
Software suite of various numerical ODE and PDE solution methods, written in Julia for the EMAT30008 (Scientific Computing) module.

## Features
Solves systems of ODEs in any number of dimensions:
![predprey](https://user-images.githubusercontent.com/45373428/165288682-bee465b2-e36f-490a-94f6-33550da8ba8a.png)

Plot Phase portraits:
![circle](https://user-images.githubusercontent.com/45373428/165285013-301a8ac7-b54f-4640-9e83-f6ccecd44daa.png)

Plot Phase portraits in 3d:
![3d-phase-portrait](https://user-images.githubusercontent.com/45373428/165290092-86cdc8bf-c770-48a3-9372-eb9c5af81927.png)


Perform numerical shooting and find limit cycles:
![predprey-phase-portrait](https://user-images.githubusercontent.com/45373428/165290265-32cf3e97-cf6e-43ee-b30a-60cf14b7ba52.png)


Find bifurcations using numerical continuation:
![predprey-continuation](https://user-images.githubusercontent.com/45373428/165290587-2d1907aa-b736-4daa-9ae0-9a6d7a624048.png)


Track parameters in ODE systems using the pseudo-arclength method:
![pseudo-arclength](https://user-images.githubusercontent.com/45373428/165290301-f1d15bb2-08b6-4034-abd3-7f8a34978895.png)



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



