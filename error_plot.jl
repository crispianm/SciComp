include("ode_solver.jl")

x0 = [1.]
t = [0, 1]
h = exp10.(range(-5,0,80))
real = exp(1)    # real value of x(1)

euler_error = Any[]
rk4_error = Any[]

for hvalue in h
    euler_sol = solve_ode(f, x0, t, euler_step, hvalue)[end][1]
    rk4_sol = solve_ode(f, x0, t, rk4_step, hvalue)[end][1]
    push!(euler_error, abs.(euler_sol.-real))
    push!(rk4_error, abs.(rk4_sol.-real))
end

using PlotlyJS


t1 = scatter(;x=h,
              y=euler_error,
              mode="lines",
              name="Euler")

t2 = scatter(;x=h,
              y=rk4_error,
              mode="lines",
              name="RK4")



layout = Layout(xaxis_type="log", xaxis_exponentformat = "power", yaxis_type="log",  yaxis_exponentformat = "power")

data = [t1, t2]

PlotlyJS.plot(data, layout)