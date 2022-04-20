include("../visualisation.jl")
include("../ode_solver.jl")
include("../numerical_shooting.jl")
include("../examples/example_functions.jl")

# function hopf2d(u, t; beta=1, sigma=-1.0)

#     u1, u2 = u
#     du1dt = beta*u1 - u2 + sigma*u1*(u1^2 + u2^2)
#     du2dt = u1 + beta*u2 + sigma*u2*(u1^2 + u2^2)

#     return [du1dt du2dt]
# end

# function hopf_sol(t; beta=1, theta=0.0)

#     u1 = √(beta) * cos.(t .+ theta)
#     u2 = √(beta) * sin.(t .+ theta)

#     return u1, u2
# end

# Plot ODE for first 40 seconds
t = 0:0.1:40
u0 = [-1 0]
plot_ode((u, t) -> hopf2d(u, t, beta=1), u0, t, ["u1" "u2"])

# Plot phase portrait
plot_phase_portrait((u, t) -> hopf2d(u, t, beta=1), u0, t, ["u1" "u2"])

# Find limit cycle
u0, T = find_limit_cycle((u, t) -> hopf2d(u, t, beta=1), [-1 0], 6.28)
println("U0: ", u0)
println("Period: ", T)


# Plot ODE analytic and numerical solutions to see if they match
t = 0:0.1:T
u1, u2 = hopf_sol(t, beta=1, theta=pi)

# numerical solution
deltat_max = 0.01
numerical_sol = solve_ode((u, t) -> hopf2d(u, t, beta=1), u0, t, "rk4", deltat_max)
x = numerical_sol[:,1]
y = numerical_sol[:,2]
x_plot = scatter(x=t, y=x, mode="lines", name="u1 - numerical", showlegend=true)
y_plot = scatter(x=t, y=y, mode="lines", name="u2 - numerical", showlegend=true)

# analytic solution
u1_plot = scatter(x=t, y=u1, mode="lines", name="u1 - analytic", showlegend=true)
u2_plot = scatter(x=t, y=u2, mode="lines", name="u2 - analytic", showlegend=true)

layout = Layout(xaxis_title = "time")

plot([x_plot, y_plot, u1_plot, u2_plot], layout)



# using Test
# function foo(x)
#     length(x)^2
# end
# println(@test foo("bar") == 9)
# println(@test foo("bar") == 10)