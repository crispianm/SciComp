include("../visualisation.jl")
include("../ode_solver.jl")
include("../numerical_shooting.jl")
include("../examples/example_functions.jl")

t = 0:0.01:100

plot_phase_portrait_3d(hopf3d, [1 0 -2], t)
# plot_phase_portrait_3d(cheng_wang, [1 0 1], t)
# plot_phase_portrait_3d(lorenz, [1 1 1], t)



# using Test
# function foo(x)
#     length(x)^2
# end
# println(@test foo("bar") == 9)
# println(@test foo("bar") == 10)