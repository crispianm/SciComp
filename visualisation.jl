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

function plot_phase_portrait()

    # Define the ODE
    # ode = @(t, u) [u[2]; -u[1] - u[2]*u[2]]

    # Define the initial conditions
    u0 = [1; 0]

    # Define the time interval
    t = linspace(0, 10, 100)

    # Plot the phase portrait
    plot_ode(ode, u0, t)

end

function  plot_nullcline()
    
        # Define the ODE
        # ode = @(t, u) [u[2]; -u[1] - u[2]*u[2]]
    
        # Define the initial conditions
        u0 = [1; 0]
    
        # Define the time interval
        t = linspace(0, 10, 100)
    
        # Plot the nullcline
        plot_ode(ode, u0, t, ["t" "u"])
    
end