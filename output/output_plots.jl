# Import Subsystems
include("../examples/example_functions.jl")
include("../finite_difference.jl")
include("../numerical_continuation.jl")
include("../numerical_shooting.jl")
include("../ode_solver.jl")
include("../visualisation.jl")
using PlotlyJS


# Set save_figures and save_figures_3d variables to true to save pngs to the output folder.
# Set save_html to true to save figures as interactive html files
save_figures = true
save_figures_3d = true
save_html = true


println("Generating output plots")

"""
Plotting 
"""

if save_figures

    ## Plot of ODE for t = 20
    t = 0:0.02:20
    u0 = [-1 0]
    hopf_plot = plot_ode(hopf2d, u0, t, ["u1 (x)" "u2 (y)"])
    savefig(hopf_plot, "./output/Hopf Plot.png")


    ## Plot 3d phase portrait
    hopf_phase_portrait = plot_phase_portrait(hopf2d, u0, t, ["u1 (x)" "u2 (y)"])
    savefig(hopf_phase_portrait, "./output/Hopf Phase Portrait.png")


    ## Plot ODE analytic and numerical solutions to see if they match
    # Find limit cycle
    u0, T = find_limit_cycle(hopf2d, [-1 0], 6.28)
    t = 0:0.01:T
    hopf_solution = hopf2d_sol(t, beta = 1, theta = pi)

    # Numerical solution
    numerical_sol = solve_ode(hopf2d, u0, t)

    x_plot = scatter(
        x = t,
        y = numerical_sol[:, 1],
        mode = "lines",
        name = "u1 (x) - numerical",
        showlegend = true,
    )
    y_plot = scatter(
        x = t,
        y = numerical_sol[:, 2],
        mode = "lines",
        name = "u2 (y) - numerical",
        showlegend = true,
    )

    # Analytic solution
    u1_plot = scatter(
        x = t,
        y = hopf_solution[:, 1],
        mode = "lines",
        name = "u1 (x)  - analytic",
        showlegend = true,
    )
    u2_plot = scatter(
        x = t,
        y = hopf_solution[:, 2],
        mode = "lines",
        name = "u2 (y)  - analytic",
        showlegend = true,
    )

    layout = Layout(xaxis_title = "time")

    sol_plot = plot([x_plot, y_plot, u1_plot, u2_plot], layout)
    savefig(sol_plot, "./output/Hopf Solution Plot.png")


    ## Finite Difference method comparison plot
    kappa, L, T, mx, mt = 1, 1, 0.5, 10, 1000

    fe_x, fe_u_j = finite_difference(u_I, kappa, L, T, mx, mt, method = "fe") # forward euler Estimate  
    be_x, be_u_j = finite_difference(u_I, kappa, L, T, mx, mt, method = "be") # backward euler Estimate
    cn_x, cn_u_j = finite_difference(u_I, kappa, L, T, mx, mt, method = "cn") # crank nicholson Estimate

    layout = Layout(xaxis_title = "x", yaxis_title = string("u(x, ", string(T), ")"))

    fd_method_comp = plot(
        [
            scatter(
                x = 0:L/250:L,
                y = u_exact(0:L/250:L, T),
                mode = "lines",
                name = "Exact Solution",
                showlegend = true,
            ),
            scatter(
                x = fe_x,
                y = fe_u_j,
                mode = "markers",
                name = "Forward Euler",
                showlegend = true,
            ),
            scatter(
                x = be_x,
                y = be_u_j,
                mode = "markers",
                name = "Backward Euler",
                showlegend = true,
            ),
            scatter(
                x = cn_x,
                y = cn_u_j,
                mode = "markers",
                name = "Crank Nicholson",
                showlegend = true,
            ),
        ],
        layout,
    )

    savefig(fd_method_comp, "./output/Finite Difference Method Comparison.png")
end



"""
3d plotting stuff

This is just to check if the visualisation and the numerical shooting methods
are working, by comparing visual output to that of known solutions. Also I
found it really fun and procrastinated a lot with it.
"""

if save_figures_3d

    ## Plot 3d phase portrait of Hopf 3d System
    t = 0:0.1:200
    u0 = [1 0 -2]
    hopf_phase_portrait_3d = plot_phase_portrait_3d(hopf3d, u0, t)
    savefig(hopf_phase_portrait_3d, "./output/3d Phase Portrait - Hopf.png")



    ## Plot 3d phase portrait of Cheng-Wang System
    t = 0:0.01:100
    u0 = [1 0 1]
    cheng_wang_phase_portrait_3d = plot_phase_portrait_3d(cheng_wang, u0, t)
    savefig(cheng_wang_phase_portrait_3d, "./output/3d Phase Portrait - Cheng-Wang.png")



    ## Plot 3d phase portrait of Lorenz System
    t = 0:0.01:100
    u0 = [1 1 1]
    lorenz_phase_portrait_3d = plot_phase_portrait_3d(lorenz, u0, t)
    savefig(lorenz_phase_portrait_3d, "./output/3d Phase Portrait - Lorenz.png")



    ## 3d Finite Difference method comparison plot
    kappa, L, T, mx, mt = 0.5, 1, 1, 10, 100

    x = range(0, stop = L, length = mx + 1)
    t = range(0.0001, stop = T, length = mt)

    u_est = ones(mt)' .* x
    v_est = t' .* ones(mx + 1)

    fe_z = Matrix{Float64}(undef, mx + 1, mt)
    be_z = Matrix{Float64}(undef, mx + 1, mt)
    cn_z = Matrix{Float64}(undef, mx + 1, mt)

    for i in range(1, mt)
        fex, fe_uj = finite_difference(u_I, kappa, L, t[i], mx, mt, method = "fe")
        fe_z[:, i] = fe_uj
        bex, be_uj = finite_difference(u_I, kappa, L, t[i], mx, mt, method = "be")
        be_z[:, i] = be_uj
        cnx, cn_uj = finite_difference(u_I, kappa, L, t[i], mx, mt, method = "cn")
        cn_z[:, i] = cn_uj
    end

    # Estimate Surfaces
    fe_estimate = scatter(
        x = vec(u_est),
        y = vec(v_est),
        z = vec(fe_z),
        name = "Forward Euler",
        mode = "markers",
        type = "scatter3d",
        marker = attr(color = "red", size = 2),
    )
    be_estimate = scatter(
        x = vec(u_est),
        y = vec(v_est),
        z = vec(be_z),
        name = "Backward Euler",
        mode = "markers",
        type = "scatter3d",
        marker = attr(color = "blue", size = 2),
    )
    cn_estimate = scatter(
        x = vec(u_est),
        y = vec(v_est),
        z = vec(cn_z),
        name = "Crank Nicholson",
        mode = "markers",
        type = "scatter3d",
        marker = attr(color = "green", size = 2),
    )


    # Exact surface
    x = range(0, stop = L, length = 100)
    t = range(0, stop = T, length = 100)

    u = ones(100)' .* x
    v = t' .* ones(100)
    z = u_exact(u, v)

    exact = surface(x = u, y = v, z = z, name = "Exact Solution", opacity = 1)

    # Plot
    layout = Layout(
        showlegend = true,
        legend = attr(x = 0, y = 1),
        legend_itemsizing = "constant",
        xaxis_title = "x",
        yaxis_title = "t",
        zaxis_title = "z",
    )

    fd_method_3d_comp = plot([exact, fe_estimate, be_estimate, cn_estimate], layout)
    savefig(fd_method_3d_comp, "./output/Finite Difference Method 3d Comparison.png")

end



"""
Save Plots as HTML files

This is so you can play with the plots and zoom in and stuff, which is pretty fun.
"""

if save_html

    savefig(hopf_plot, "./output/Hopf Plot.html")
    savefig(hopf_phase_portrait, "./output/Hopf Phase Portrait.html")
    savefig(sol_plot, "./output/Hopf Solution Plot.html")
    savefig(fd_method_comp, "./output/Finite Difference Method Comparison.html")
    savefig(hopf_phase_portrait_3d, "./output/3d Phase Portrait - Hopf.html")
    savefig(cheng_wang_phase_portrait_3d, "./output/3d Phase Portrait - Cheng-Wang.html")
    savefig(lorenz_phase_portrait_3d, "./output/3d Phase Portrait - Lorenz.html")
    savefig(fd_method_3d_comp, "./output/Finite Difference Method 3d Comparison.html")
end

println("Output generated.")