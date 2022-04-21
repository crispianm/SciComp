using PlotlyJS
using Test
include("../ode_solver.jl")
include("../numerical_shooting.jl")
include("../examples/example_functions.jl")
include("../visualisation.jl")

# Plot of ODE for t = 40
t = 0:0.1:40
u0 = [-1 0]
hopf_plot = plot_ode((u, t) -> hopf2d(u, t, beta=1), u0, t, ["x" "y"])
savefig(hopf_plot, "./output/hopf_plot.png")

# Plot phase portrait
hopf_phase_portrait = plot_phase_portrait((u, t) -> hopf2d(u, t, beta=1), u0, t, ["u1" "u2"])
savefig(hopf_phase_portrait, "./output/hopf_phase_portrait.png")

# Find limit cycle
u0, T = find_limit_cycle((u, t) -> hopf2d(u, t, beta=1), [-1 0], 6.28)
println("u0: ", u0)
println("Period: ", T)


# Plot ODE analytic and numerical solutions to see if they match
t = 0:0.1:T
u1, u2 = hopf2d_sol(t, beta=1, theta=pi)

# Numerical solution
deltat_max = 0.01
numerical_sol = solve_ode((u, t) -> hopf2d(u, t, beta=1), u0, t, "rk4", deltat_max)
x = numerical_sol[:,1]
y = numerical_sol[:,2]
x_plot = scatter(x=t, y=x, mode="lines", name="u1 - numerical", showlegend=true)
y_plot = scatter(x=t, y=y, mode="lines", name="u2 - numerical", showlegend=true)

# Analytic solution
u1_plot = scatter(x=t, y=u1, mode="lines", name="u1 - analytic", showlegend=true)
u2_plot = scatter(x=t, y=u2, mode="lines", name="u2 - analytic", showlegend=true)

layout = Layout(xaxis_title = "time")

sol_plot = plot([x_plot, y_plot, u1_plot, u2_plot], layout)
savefig(sol_plot, "./output/sol_plot.png")





x0 = [1]
t = [0, 1]

estimate = solve_ode(f, x0, t, "rk4", 0.001)[end][1]

println(@test isapprox(estimate, ℯ))
println(@test isapprox(estimate, ℯ))
println(@test isapprox(estimate, ℯ))
println(@test isapprox(estimate, ℯ))
println(@test isapprox(estimate, ℯ))
println(@test isapprox(estimate, ℯ))
println(@test isapprox(estimate, ℯ))
println(@test isapprox(estimate, ℯ))
println(@test isapprox(estimate, ℯ))
println(@test isapprox(estimate, ℯ))
println(@test isapprox(estimate, ℯ))