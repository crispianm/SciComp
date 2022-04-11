using PlotlyJS
include("ode_solver.jl")

function plot_ode(ode, u0, t, labels=["t" "x"])

    deltat_max = 0.01
    solution = solve_ode(ode, u0, t, rk4_step, deltat_max)
    x = solution[:,1]
    y = solution[:,2]

    # Create traces
    rk4_x = scatter(x=t, y=x, mode="lines", name=labels[1], showlegend=true)
    rk4_y = scatter(x=t, y=y, mode="lines", name=labels[2], showlegend=true)


    layout = Layout(
        xaxis_title = "time",
        yaxis_title = "x")

    plot([rk4_x, rk4_y], layout)

end


function plot_ode_3d(ode, u0, t, labels=["u1" "u2" "u3"])

    deltat_max = 0.01
    solution = solve_ode(ode, u0, t, rk4_step, deltat_max)
    x = solution[:,1]
    y = solution[:,2]
    z = solution[:,3]

    # Create traces
    rk4_x = scatter(x=t, y=x, mode="lines", name=labels[1], showlegend=true)
    rk4_y = scatter(x=t, y=y, mode="lines", name=labels[2], showlegend=true)
    rk4_z = scatter(x=t, y=z, mode="lines", name=labels[3], showlegend=true)


    layout = Layout(
        xaxis_title = "time",
        yaxis_title = "x")

    plot([rk4_x, rk4_y, rk4_z], layout)

end


function plot_phase_portrait(ode, u0, t, labels=["x" "y"])

    deltat_max = 0.01
    solution = solve_ode(ode, u0, t, rk4_step, deltat_max)
    x = solution[:,1]
    y = solution[:,2]

    # Create traces
    rk4 = scatter(x=x, y=y, mode="lines", name="rk4 approximation of ODE", showlegend=true)

    layout = Layout(
        xaxis_title = labels[1],
        yaxis_title = labels[2])

    plot([rk4], layout)

end


function plot_phase_portrait_3d(ode, u0, t, labels=["u1" "u2" "u3"])

    deltat_max = 0.01
    solution = solve_ode(ode, u0, t, rk4_step, deltat_max)
    x = solution[:,1]
    y = solution[:,2]
    z = solution[:,3]
    
    # Create traces
    # rk4 = scatter_3d(x=x, y=y, z=z, mode="lines", name="rk4 approximation of ODE", showlegend=true)

    # layout = Layout()

    # plot([rk4], layout)

    plot(
        scatter(
            x=x,
            y=y,
            z=z,
            mode="lines",
            type="scatter3d"
        ),
        Layout(
            margin=attr(l=0, r=0, b=0, t=0),
            xaxis_title = labels[1],
            yaxis_title = labels[2],
            zaxis_title = labels[3])
        )

end


function  plot_nullcline()
    
    whirrr
    
end