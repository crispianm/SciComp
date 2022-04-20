using PlotlyJS
include("ode_solver.jl")

function plot_ode(ode, u0, t, labels=["t" "x"], arg...)

    """
    Plots a provided ODE, ode, along time input t with initia condition(s) u0.
    
    Parameters:
        ode (function): Function which returns a singular value or 1 x n matrix of values.
        u0 (matrix): Matrix of initial values in the 1 x n form, eg: [1] or [1 1].
        t (array or range): Time values to solve between. The first element of t should be the initial value, 0.
        labels (array, optional): Labels for the x and y axes.
        labels (array, optional): Labels for the x and y components to plot.
        arg (list, optional): Arguments to pass to f.
    
    Returns:
        A plot of the solutions for x and y at each time value in t
    """

    deltat_max = 0.01
    solution = solve_ode(ode, u0, t, "rk4", deltat_max)
    x = solution[:,1]
    y = solution[:,2]

    # Create traces
    rk4_x = scatter(x=t, y=x, mode="lines", name=labels[1], showlegend=true)
    rk4_y = scatter(x=t, y=y, mode="lines", name=labels[2], showlegend=true)


    layout = Layout(
        xaxis_title = "time",
        width=700, height=350,
        )

    plot([rk4_x, rk4_y], layout)

end


function plot_ode_3d(ode, u0, t, labels=["u1" "u2" "u3"], arg...)

    """
    Plots a provided ODE, ode, in 3d along time input t with initia condition(s) u0.
    
    Parameters:
        ode (function): Function which returns a singular value or 1 x n matrix of values.
        u0 (matrix): Matrix of initial values in the 1 x n form, eg: [1] or [1 1].
        t (array or range): Time values to solve between. The first element of t should be the initial value, 0.
        labels (array, optional): Labels for the x, y, and z components to plot.
        arg (list, optional): Arguments to pass to f.
    
    Returns:
        x_series: List of solutions for x at each time value in t
    """

    deltat_max = 0.01
    solution = solve_ode(ode, u0, t, "rk4", deltat_max)
    x = solution[:,1]
    y = solution[:,2]
    z = solution[:,3]

    # Create traces
    rk4_x = scatter(x=t, y=x, mode="lines", name=labels[1], showlegend=true)
    rk4_y = scatter(x=t, y=y, mode="lines", name=labels[2], showlegend=true)
    rk4_z = scatter(x=t, y=z, mode="lines", name=labels[3], showlegend=true)


    layout = Layout(
        xaxis_title = "time",
        width=700, height=350,
        )

    plot([rk4_x, rk4_y, rk4_z], layout)

end


function plot_phase_portrait(ode, u0, t, labels=["x" "y"], arg...)

    """
    Plots a provided ODE, ode, along time input t with initia condition(s) u0.
    
    Parameters:
        ode (function): Function which returns a singular value or 1 x n matrix of values.
        u0 (matrix): Matrix of initial values in the 1 x n form, eg: [1] or [1 1].
        t (array or range): Time values to solve between. The first element of t should be the initial value, 0.
        labels (array, optional): Labels for the x and y axes.
        arg (list, optional): Arguments to pass to f.
    
    Returns:
        x_series: List of solutions for x at each time value in t
    """
    
    deltat_max = 0.01
    solution = solve_ode(ode, u0, t, "rk4", deltat_max)
    x = solution[:,1]
    y = solution[:,2]

    # Create trace
    rk4 = scatter(
        x=x,
        y=y,
        mode="lines",
        name="rk4 approximation of ODE",
        showlegend=true
        )

    layout = Layout(
        xaxis_title = labels[1],
        yaxis_title = labels[2],
        width=700, height=350,
        )

    plot([rk4], layout)

end


function plot_phase_portrait_3d(ode, u0, t, labels=["u1" "u2" "u3"], arg...)

    deltat_max = 0.01
    solution = solve_ode(ode, u0, t, "rk4", deltat_max)
    x = solution[:,1]
    y = solution[:,2]
    z = solution[:,3]

    # Create trace
    rk4 = scatter(
        x=x,
        y=y,
        z=z,
        mode="lines",
        type="scatter3d"
        )

    layout = Layout(
        margin=attr(l=0, r=0, b=0, t=0),
        xaxis_title = labels[1],
        yaxis_title = labels[2],
        zaxis_title = labels[3],
        width=700, height=350,
        )

    plot([rk4], layout)

end

