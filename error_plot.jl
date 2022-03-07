include("ode_solver.jl")

t = 0:1:3
# t = 1:1:3

solve_ode(f, 1, t, euler_step, 0.001)