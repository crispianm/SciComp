using PlotlyJS
include("../ode_solver.jl")

function f(x,t)
    return x
end

x0 = [1]
t = [0 1]
Δt_max = exp10.(range(-5,0,80))
real = ℯ

println("For Δt = 0.001: ")
# Euler estimate of x(1)
solution = solve_ode(f, x0, t; method="euler", Δt=0.001)
println("\tEuler approximation = \t\t", solution[end][1])

# Ralston4 estimate of x(1)
solution = solve_ode(f, x0, t; method="ralston4", Δt=0.001)
println("\tRalston4 approximation = \t", solution[end][1])

# Heun3 estimate of x(1)
solution = solve_ode(f, x0, t; method="heun3", Δt=0.001)
println("\tHeun3 approximation = \t\t", solution[end][1])

# 3/8 Rule estimate of x(1)
solution = solve_ode(f, x0, t; method="three_eighths_rule", Δt=0.001)
println("\t3/8 Rule approximation = \t", solution[end][1])

# RK4 estimate of x(1)
solution = solve_ode(f, x0, t; method="rk4", Δt=0.001)
println("\tRK4 approximation = \t\t", solution[end][1])

euler_error = Any[]
ralston4_error = Any[]
heun3_error = Any[]
three_eighths_rule_error = Any[]
rk4_error = Any[]

for value in Δt_max
    euler_sol = solve_ode(f, x0, t; method="euler", Δt=value)[end][1];
    ralston4_sol = solve_ode(f, x0, t; method="ralston4", Δt=value)[end][1];
    heun3_sol = solve_ode(f, x0, t; method="heun3", Δt=value)[end][1];
    three_eighths_rule_sol = solve_ode(f, x0, t; method="three_eighths_rule", Δt=value)[end][1];
    rk4_sol = solve_ode(f, x0, t; method="rk4", Δt=value)[end][1];

    push!(euler_error, abs.(euler_sol .- real))
    push!(ralston4_error, abs.(ralston4_sol .- real))
    push!(heun3_error, abs.(heun3_sol .- real))
    push!(three_eighths_rule_error, abs.(three_eighths_rule_sol .- real))   
    push!(rk4_error, abs.(rk4_sol .- real))
end


t1 = scatter(x=Δt_max,
              y=euler_error,
              mode="lines",
              name="Euler")
              
t2 = scatter(x=Δt_max,
            y=ralston4_error,
            mode="lines",
            name="Ralston 4")

t3 = scatter(x=Δt_max,
              y=heun3_error,
              mode="lines",
              name="Heun 3")


t4 = scatter(x=Δt_max,
              y=three_eighths_rule_error,
              mode="lines",
              name="3/8 rule")

t5 = scatter(x=Δt_max,
              y=rk4_error,
              mode="lines",
              name="RK4")

layout = Layout(
    xaxis_type="log",
    xaxis_exponentformat="power",
    xaxis_title="Δt",
    yaxis_type="log",
    yaxis_exponentformat="power",
    yaxis_title="error",
    width=700, height=350,
)
data = [t1, t2, t3, t4, t5]

plot(data, layout)