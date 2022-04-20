include("../visualisation.jl")
include("../ode_solver.jl")
include("../numerical_shooting.jl")
include("../examples/example_functions.jl")

"""
This file is just to check if the visualisation and the numerical shooting are working,
by visually comparing output to that of known solutions.
"""


# Hopf System

t = 0:0.1:200
u0 = [1 0 -2]
hopf_phase_portrait_3d = plot_phase_portrait_3d(hopf3d, u0, t)
savefig(hopf_phase_portrait_3d, "./output/hopf_phase_portrait_3d.png")

# Cheng-Wang System

t = 0:0.01:100
u0 = [1 0 1]
cheng_wang_phase_portrait_3d = plot_phase_portrait_3d(cheng_wang, u0, t)
savefig(cheng_wang_phase_portrait_3d, "./output/cheng_wang_phase_portrait_3d.png")

# Lorenz System

t = 0:0.01:100
u0 = [1 1 1]
lorenz_phase_portrait_3d = plot_phase_portrait_3d(lorenz, u0, t)
savefig(lorenz_phase_portrait_3d, "./output/lorenz_phase_portrait_3d.png")
