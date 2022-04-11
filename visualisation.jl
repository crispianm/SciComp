using PlotlyJS
include("ode_solver.jl")

function plot_ode(ode, u0, t, labels=["t" "x"])

    deltat_max = 0.01
    solution = solve_ode(ode, u0, t, rk4_step, deltat_max)[:,1]

    # Create traces
    rk4 = scatter(x=t, y=solution, mode="lines", name="rk4 approximation of ODE", showlegend=true)

    layout = Layout(
        xaxis_title = labels[1],
        yaxis_title = labels[2])

    plot([rk4], layout)

end

function plot_phase_portrait(ode, u0, t)

    deltat_max = 0.01
    solution = solve_ode(ode, u0, t, rk4_step, deltat_max)
    x = solution[:,1]
    y = solution[:,2]

    # Create traces
    rk4 = scatter(x=x, y=y, mode="lines", name="rk4 approximation of ODE", showlegend=true)

    layout = Layout()

    plot([rk4], layout)

end

function  plot_nullcline()
    
    whirrr
    
end