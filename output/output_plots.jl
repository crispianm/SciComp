# Import All Subsystems
include("../examples/example_functions.jl")
include("../src/FiniteDifference.jl")
include("../src/NumericalContinuation.jl")
include("../src/NumericalShooting.jl")
include("../src/ODESolver.jl")
include("../src/Visualisation.jl")
using PlotlyJS


"""
Creates plotly graphs of various ODEs and PDEs from the example functions
file, showing some of the capabilities of the software. The priority was on 
making the plots as clear as possible, occasionally at some code readability
or repetition cost (see Continuation Plot).

It takes about 40-50 seconds to run if all figures are chosen to be generated.

Set save_figures and save_figures_3d variables to true to save pngs to the output folder.
Set save_html to true to save figures there as interactive html files, which are loads of fun.

"""


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
    hopf_plot = plot_ode(hopf2d, u0, t, labels=["u1 (x)" "u2 (y)"])
    savefig(hopf_plot, "./output/Hopf Plot.png")


    ## Plot Hopf phase portrait
    hopf_phase_portrait = plot_phase_portrait(hopf2d, u0, t, labels=["u1 (x)" "u2 (y)"])
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


    ## Continuation Plot
    # Excercise 1
    algebraic_pv_1, algebraic_c_1 = continuation(
        algebraic,
        [1 1],
        6,
        "c",
        -2:0.05:2,
        method = "natural_parameter_continuation",
        discretisation = "none",
    )
    hopf2d_pv_1, hopf2d_c_1 = continuation(
        hopf2d,
        [-1 0],
        6,
        "beta",
        2:-0.05:0,
        method = "natural_parameter_continuation",
        discretisation = "shooting",
    )
    hopf2d_mod_pv_1, hopf2d_mod_c_1 = continuation(
        hopf2d_modified,
        [-1 0],
        6,
        "beta",
        2:-0.05:-1,
        method = "natural_parameter_continuation",
        discretisation = "shooting",
    )

    # Excercise 2
    algebraic_pv_2, algebraic_c_2 = continuation(
        algebraic,
        [1 1],
        6,
        "c",
        -2:0.05:2,
        method = "pseudo_arclength",
        discretisation = "none",
    )
    hopf2d_pv_2, hopf2d_c_2 = continuation(
        hopf2d,
        [-1 0],
        6,
        "beta",
        2:-0.05:0,
        method = "pseudo_arclength",
        discretisation = "shooting",
    )
    hopf2d_mod_pv_2, hopf2d_mod_c_2 = continuation(
        hopf2d_modified,
        [-1 0],
        6,
        "beta",
        2:-0.05:-1,
        method = "pseudo_arclength",
        discretisation = "shooting",
    )

    # Create traces
    p1_u1 = scatter(
        x = algebraic_pv_1,
        y = algebraic_c_1[:, 1],
        mode = "lines",
        name = "u1",
        showlegend = true,
        line = attr(color = "blue"),
    )
    p1_u2 = scatter(
        x = algebraic_pv_1,
        y = algebraic_c_1[:, 2],
        mode = "lines",
        name = "u2",
        showlegend = true,
        line = attr(color = "red"),
    )

    p2_u1 = scatter(
        x = hopf2d_pv_1,
        y = hopf2d_c_1[:, 1],
        mode = "lines",
        name = "u1",
        showlegend = false,
        line = attr(color = "blue"),
    )
    p2_u2 = scatter(
        x = hopf2d_pv_1,
        y = hopf2d_c_1[:, 2],
        mode = "lines",
        name = "u2",
        showlegend = false,
        line = attr(color = "red"),
    )

    p3_u1 = scatter(
        x = hopf2d_mod_pv_1,
        y = hopf2d_mod_c_1[:, 1],
        mode = "lines",
        name = "u1",
        showlegend = false,
        line = attr(color = "blue"),
    )
    p3_u2 = scatter(
        x = hopf2d_mod_pv_1,
        y = hopf2d_mod_c_1[:, 2],
        mode = "lines",
        name = "u2",
        showlegend = false,
        line = attr(color = "red"),
    )

    p4_u1 = scatter(
        x = algebraic_pv_2,
        y = algebraic_c_2[:, 1],
        mode = "lines",
        name = "u1",
        showlegend = false,
        line = attr(color = "blue"),
    )
    p4_u2 = scatter(
        x = algebraic_pv_2,
        y = algebraic_c_2[:, 2],
        mode = "lines",
        name = "u2",
        showlegend = false,
        line = attr(color = "red"),
    )

    p5_u1 = scatter(
        x = hopf2d_pv_2,
        y = hopf2d_c_2[:, 1],
        mode = "lines",
        name = "u1",
        showlegend = false,
        line = attr(color = "blue"),
    )
    p5_u2 = scatter(
        x = hopf2d_pv_2,
        y = hopf2d_c_2[:, 2],
        mode = "lines",
        name = "u2",
        showlegend = false,
        line = attr(color = "red"),
    )

    p6_u1 = scatter(
        x = hopf2d_mod_pv_2,
        y = hopf2d_mod_c_2[:, 1],
        mode = "lines",
        name = "u1",
        showlegend = false,
        line = attr(color = "blue"),
    )
    p6_u2 = scatter(
        x = hopf2d_mod_pv_2,
        y = hopf2d_mod_c_2[:, 2],
        mode = "lines",
        name = "u2",
        showlegend = false,
        line = attr(color = "red"),
    )


    # Plot
    p1 = plot(
        [p1_u1, p1_u2],
        Layout(title = "Algebraic NPC", xaxis_title = "c", yaxis_title = ""),
    )
    p2 = plot(
        [p2_u1, p2_u2],
        Layout(title = "Hopf NPC", xaxis_title = "beta", yaxis_title = ""),
    )
    p3 = plot(
        [p3_u1, p3_u2],
        Layout(title = "Modified Hopf NPC", xaxis_title = "beta", yaxis_title = ""),
    )
    p4 = plot(
        [p4_u1, p4_u2],
        Layout(title = "Algebraic PA", xaxis_title = "c", yaxis_title = ""),
    )
    p5 = plot(
        [p5_u1, p5_u2],
        Layout(title = "Hopf PA", xaxis_title = "beta", yaxis_title = ""),
    )
    p6 = plot(
        [p6_u1, p6_u2],
        Layout(title = "Modified Hopf PA", xaxis_title = "beta", yaxis_title = ""),
    )


    ncm_comp_plot = [p1 p2 p3; p4 p5 p6]
    relayout!(
        ncm_comp_plot,
        height = 500,
        width = 800,
        title_text = "Numerical Continuation Method Comparison",
    )
    savefig(ncm_comp_plot, "./output/Continuation Method Comparison.png")
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

    x = range(0, stop=L, length=mx+1)
    t = range(0, stop=T, length=mt) .+ 0.0001
    
    u_est = ones(mt)' .* x
    v_est = t' .* ones(mx+1)
    
    fe_z = Matrix{Float64}(undef, mx+1, mt)
    be_z = Matrix{Float64}(undef, mx+1, mt)
    cn_z = Matrix{Float64}(undef, mx+1, mt)
    
    for i in 1:mt
        fex, fe_uj = finite_difference(u_I, kappa, L, t[i], mx, mt, method="fe")
        fe_z[:,i] = fe_uj
        bex, be_uj = finite_difference(u_I, kappa, L, t[i], mx, mt, method="be")
        be_z[:,i] = be_uj
        cnx, cn_uj = finite_difference(u_I, kappa, L, t[i], mx, mt, method="cn")
        cn_z[:,i] = cn_uj
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
    savefig(ncm_comp_plot, "./output/Continuation Method Comparison.html")
    savefig(hopf_phase_portrait_3d, "./output/3d Phase Portrait - Hopf.html")
    savefig(cheng_wang_phase_portrait_3d, "./output/3d Phase Portrait - Cheng-Wang.html")
    savefig(lorenz_phase_portrait_3d, "./output/3d Phase Portrait - Lorenz.html")
    savefig(fd_method_3d_comp, "./output/Finite Difference Method 3d Comparison.html")
end

println("Output generated.")
