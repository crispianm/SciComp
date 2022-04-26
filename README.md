## About
Software suite of various numerical ODE and PDE solution methods, written in Julia for the Engineering Mathematics Scientific Computing module.

## Features
**Solves systems of ODEs in any number of dimensions:**
![predprey](https://user-images.githubusercontent.com/45373428/165288682-bee465b2-e36f-490a-94f6-33550da8ba8a.png)

**Plot Phase portraits:**
![circle](https://user-images.githubusercontent.com/45373428/165285013-301a8ac7-b54f-4640-9e83-f6ccecd44daa.png)

**Plot Phase portraits in 3d:**
![3d-phase-portrait](https://user-images.githubusercontent.com/45373428/165290092-86cdc8bf-c770-48a3-9372-eb9c5af81927.png)


**Perform numerical shooting and find limit cycles:**
![predprey-phase-portrait](https://user-images.githubusercontent.com/45373428/165292376-9de1c721-2716-47f2-b2a1-4bbe93de6213.png)


**Find bifurcations using numerical continuation:**
![predprey-continuation](https://user-images.githubusercontent.com/45373428/165292421-1a154555-78b6-4745-97f5-1e2f21dd5da0.png)


**Track parameters in ODE systems using the pseudo-arclength method:**
![pseudo-arclength](https://user-images.githubusercontent.com/45373428/165290301-f1d15bb2-08b6-4034-abd3-7f8a34978895.png)

**Analyse the effects of the diffusion PDE using finite difference:**
![fd](https://user-images.githubusercontent.com/45373428/165399233-78cdb2fa-37c5-43e9-b03a-bad8be84be21.png)


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
Note: This does not work in IJulia notebooks.

**Simple Example**
```julia
using .SciCompToolkit

function ode(u,t)
    u1, u2 = u
    du1dt = u2.^2 - u1.^2
    du2dt = -u1
    return [du1dt du2dt]
end

plot_ode(ode, [1 1], 0:0.1:5)
```

