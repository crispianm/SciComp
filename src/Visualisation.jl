include("./FiniteDifference.jl")
include("./NumericalContinuation.jl")
include("./NumericalShooting.jl")
include("./ODESolver.jl")
using PlotlyJS


function plot_ode(ode, u0, t; labels = [string("u",i) for i in 1:length(u0)], arg...)

    """
    Plots a provided n-d ODE, ode, along time input t with initial condition(s) u0.
    
    Parameters:
        ode (function): Function which returns a singular value or 1 x n matrix of values.
        u0 (matrix): Matrix of initial values in the 1 x n form, eg: [1] or [1 1].
        t (array or range): Time values to solve between. The first element of t should be the initial value, 0.
        labels (array, optional): Labels for the x and y components to plot.
        arg (list, optional): Arguments to pass to f.
    
    Returns:
        A plot of the solutions for x and y at each time value in t.
    """

    if size(u0)[2] != length(labels)
        error("Please enter the same number of labels as initial conditions.")
    end

    solution = solve_ode(ode, u0, t; arg...)

    plots = GenericTrace[]
    for i in 1:size(u0)[2]

        # Calculate solution
        y = solution[:, i]

        # Create trace and append
        conditions = scatter(x = t, y = y, mode = "lines", name = labels[i], showlegend = true)
        push!(plots, conditions)

    end

    layout = Layout(
        xaxis_title = "time",
        title = string("Time series plot of ODE system: ", ode);
        arg...
    )

    plot(plots, layout)

end


function plot_phase_portrait(ode, u0, t; labels = ["u1" "u2"], arg...)

    """
    Plots a 2d phase portrait for a provided ODE, ode, along time input t
    with initia condition(s) u0.
    
    Parameters:
        ode (function): Function which returns a singular value or 1 x n matrix of values.
        u0 (matrix): Matrix of initial values in the 1 x n form, eg: [1] or [1 1].
        t (array or range): Time values to solve between. The first element of t should be the initial value, 0.
        labels (array, optional): Labels for the x and y axes.
        arg (list, optional): Arguments to pass to f.
    
    Returns:
        A 2d plot of the ode's phase portrait.
    """

    solution = solve_ode(ode, u0, t; arg...)

    # Create trace
    conditions = scatter(
        x = solution[:, 1],
        y = solution[:, 2],
        mode = "lines",
    )

    layout = Layout(
        xaxis_title = labels[1],
        yaxis_title = labels[2],
        title = string("Phase portrait of ODE system: ", ode);
        arg...
    )

    plot([conditions], layout)

end


function plot_phase_portrait_3d(ode, u0, t; labels = ["u1" "u2" "u3"], arg...)

    """
    Plots a 3d phase portrait for a provided ODE, ode, along time input t
    with initia condition(s) u0.
    
    Parameters:
        ode (function): Function which returns a singular value or 1 x n matrix of values.
        u0 (matrix): Matrix of initial values in the 1 x n form, eg: [1] or [1 1].
        t (array or range): Time values to solve between. The first element of t should be the initial value, 0.
        labels (array, optional): Labels for the x and y axes.
        arg (list, optional): Arguments to pass to f.
    
    Returns:
        A 3d plot of the ode's phase portrait.
    """

    solution = solve_ode(ode, u0, t; arg...)

    # Create trace
    conditions = scatter(
        x = solution[:, 1],
        y = solution[:, 2],
        z = solution[:, 3],
        mode = "lines",
        type = "scatter3d",
    )

    layout = Layout(
        margin = attr(l = 0, r = 0, b = 0),
        xaxis_title = labels[1],
        yaxis_title = labels[2],
        zaxis_title = labels[3],
        title = string("3d phase portrait of ODE system: ", ode);
        arg...
    )

    plot([conditions], layout)

end


function plot_continuation(ode, u0, T, parameter, par_values; labels = ["u1" "u2"], arg...)

    """
    Plots a provided ODE, ode, along parameter input 'parameter' with initia condition(s) u0.
    
    Parameters:
            ode (function): Function which returns a singular value or 1 x n matrix of values.
            u0 (matrix): Matrix of initial values in the 1 x n form, eg: [1] or [1 1].
            T (float): Initial guess for the period.
            parameter (string): The parameter in the system to vary.
                Allowable parameter inputs:
                    "a", "alpha", "b", "beta", "c", "d", "delta", "sigma"
            par_values (range): Parameter values to solve between, made with colons, eg: 0:0.1:2.
            arg (list, optional): Arguments to pass to f or continuation e.g method or discretisation.
    
    Returns:
        A plot of the initial conditions u1 and u2 at parameter value in par_values.
    """

    par_values, conditions = continuation(ode, u0, T, parameter, par_values; arg...)

    # Create traces
    u1 = scatter(x = par_values, y = conditions[:, 1], mode = "lines", name = labels[1])
    u2 = scatter(x = par_values, y = conditions[:, 2], mode = "lines", name = labels[2])


    layout = Layout(
        xaxis_title = string("Parameter: ", parameter),
        title = string("Numerical Continuation of ODE system: ", ode);
        arg...
    )

    plot([u1, u2], layout)

end

