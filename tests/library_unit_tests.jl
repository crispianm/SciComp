using .SciCompToolkit
using PlotlyJS
using Test
using ProgressMeter


"""
    Test script to run the ODE solver for a variety of differential equations
    and check them against explicit solutions.

"""


#  Set save_figures and save_figures_3d variables to true to save pngs to the output folder.
save_figures = true
save_figures_3d = true


"""
Run the tests !
"""

test = ProgressUnknown("Performing System Unit Tests: ", spinner = true)

@testset verbose = true "System Unit Tests" begin

    @testset verbose = true "ODE Solver" begin

        @testset verbose = true "Input Tests" begin

            # test error is thrown if incorrect t is input
            @test_throws ErrorException solve_ode(f2, [1], 1:0.1:2, method = "rk4")
            @test_throws ErrorException solve_ode(f2, [1], [1], method = "rk4")
            @test_throws ErrorException solve_ode(f2, [1], [1 2; 3 4], method = "rk4")


            t = 0:0.1:1
            # test error is thrown if x0 is not a matrix
            @test_throws ErrorException solve_ode(f2, 1, t, method = "rk4")
            @test_throws ErrorException solve_ode(f2, [1; 0], t, method = "rk4")
            @test_throws ErrorException solve_ode(f2, [1, 0], t, method = "rk4")


            # test error is thrown if x0 is not a 1xn matrix
            @test_throws ErrorException solve_ode(f2, [1 2; 3 4], t, method = "rk4")


            # test error is thrown if x0 is not correct length
            @test_throws ErrorException solve_ode(f2, [1], t, method = "rk4")

            ProgressMeter.next!(test)

        end

        @testset verbose = true "Output Tests" begin

            # test if solve_ode estimates a simple ODE correctly
            x0 = [1]
            t = 0:1
            e_estimate = solve_ode(f, x0, t, method = "rk4")[end][1]
            @test isapprox(e_estimate, ℯ)
            ProgressMeter.next!(test)


            # test if solve_ode estimates systems of ODEs correctly
            x0 = [1 0]
            t = 0:0.1:1
            f2_sol = f2_solution(x0, t)
            f2_numerical_sol = solve_ode(f2, x0, t, method = "rk4")
            @test all(isapprox.(f2_numerical_sol, f2_sol, atol = 1e-6))
            ProgressMeter.next!(test)

        end
    end

    @testset verbose = true "Numerical Shooting" begin

        @testset verbose = true "Input Tests" begin

            # test error is thrown if T is not an integer or float
            T = 10
            bad_T = 1:0.1:2
            @test_throws ErrorException find_limit_cycle(f2, [1], bad_T)


            # test error is thrown if u0 is not a matrix
            @test_throws ErrorException find_limit_cycle(f2, 1, T)
            @test_throws ErrorException find_limit_cycle(f2, [1; 0], T)
            @test_throws ErrorException find_limit_cycle(f2, [1, 0], T)


            # test error is thrown if u0 is not a 1xn matrix
            @test_throws ErrorException find_limit_cycle(f2, [1 2; 3 4], T)


            # test error is thrown if u0 is not correct length
            @test_throws ErrorException find_limit_cycle(f2, [1], 10)


            # test error is thrown if phase_index is a positive integer
            @test_throws ErrorException find_limit_cycle(f2, [1 1], 10, phase_index = -1)
            @test_throws ErrorException find_limit_cycle(f2, [1 1], 10, phase_index = [3])
            @test_throws ErrorException find_limit_cycle(f2, [1 1], 10, phase_index = 2.5)
            @test_throws ErrorException find_limit_cycle(f2, [1 1], 10, phase_index = "0")

            ProgressMeter.next!(test)

        end

        @testset verbose = true "Output Tests" begin

            # test if found limit cycle matches the analytical solution
            u0, T = find_limit_cycle(hopf2d, [-1 0], 6)
            @test isapprox(T, 2 * pi)
            ProgressMeter.next!(test)


            # test if solve_ode estimates a Hopf ODE correctly
            t = 0:0.1:T
            hopf_solution = hopf2d_sol(t, beta = 1, theta = pi) # adjusted for phase
            hopf_numerical_sol = solve_ode(hopf2d, u0, t, method = "rk4")
            @test all(isapprox.(hopf_numerical_sol, hopf_solution, atol = 1e-6))
            ProgressMeter.next!(test)


            # test if arguments are being passed to the function correctly
            u0_0, T = find_limit_cycle(hopf2d, [-1 0], 6)
            u0_1, T = find_limit_cycle(hopf2d, [-1 0], 6, beta = 2)
            u0_3, T = find_limit_cycle(hopf2d, [-1 0], 6, beta = 2, sigma = -1.2)
            u0_2, T = find_limit_cycle(hopf2d, [-1 0], 6, sigma = -1.2)
            @test u0_0 != u0_1
            @test u0_0 != u0_2
            @test u0_0 != u0_3

            ProgressMeter.next!(test)

        end

        @testset verbose = true "Higher Dimensional Output" begin

            # test if found limit cycle matches the analytical solution
            u0, T = find_limit_cycle(hopf3d, [1 1 1], 6)
            @test isapprox(T, 2 * pi)
            ProgressMeter.next!(test)


            u0, T = find_limit_cycle(hopf3d, [1 1 1], 10)
            hopf_solution = hopf3d_sol(u0, 0:0.1:T; theta = 0) # adjusted for phase
            hopf_numerical_sol = solve_ode(hopf3d, u0, 0:0.1:T, method = "rk4")
            @test all(isapprox.(hopf_numerical_sol, hopf_solution, atol = 1e-5))
            ProgressMeter.next!(test)

        end
    end

    @testset verbose = true "Numerical Continuation" begin

        @testset verbose = true "Input Tests" begin

            # test error is thrown if T is not an integer or float
            bad_T = 1:0.1:2
            @test_throws ErrorException continuation(hopf2d, [1 1], bad_T, "beta", 0:0.01:2)


            T = 10
            # test error is thrown if u0 is not a matrix
            @test_throws ErrorException continuation(hopf2d, 1, 6, "beta", 0:0.01:2)
            @test_throws ErrorException continuation(hopf2d, [1; 0], 6, "beta", 0:0.01:2)
            @test_throws ErrorException continuation(hopf2d, [1, 0], 6, "beta", 0:0.01:2)


            # test error is thrown if u0 is not a 1xn matrix
            @test_throws ErrorException continuation(
                hopf2d,
                [1 2; 3 4],
                6,
                "beta",
                0:0.01:2,
            )


            # test error is thrown if u0 is not correct length
            @test_throws ErrorException continuation(hopf2d, [1], 6, "beta", 0:0.01:2)


            # test error is thrown if method input is incorrect
            @test_throws ErrorException continuation(
                hopf2d,
                [1 1],
                6,
                "beta",
                0:0.01:2,
                method = 1,
            )
            @test_throws ErrorException continuation(
                hopf2d,
                [1 1],
                6,
                "beta",
                0:0.01:2,
                method = "pseudo_arclengths",
            )


            # test error is thrown if discretisation input is incorrect
            @test_throws ErrorException continuation(
                hopf2d,
                [1 1],
                6,
                "beta",
                0:0.01:2,
                discretisation = 2,
            )
            @test_throws ErrorException continuation(
                hopf2d,
                [1 1],
                6,
                "beta",
                0:0.01:2,
                discretisation = "shoot",
            )


            # test error is thrown if par_values input is incorrect
            @test_throws ErrorException continuation(hopf2d, [1 1], 6, "beta", [0 2])
            @test_throws ErrorException continuation(
                hopf2d,
                [1 1],
                6,
                "beta",
                [0.1 0.5 1.0 1.5 2.0],
            )
            @test_throws ErrorException continuation(hopf2d, [1 1], 6, "beta", 5)
            @test_throws ErrorException continuation(hopf2d, [1 1], 6, "beta", [5])
            @test_throws ErrorException continuation(hopf2d, [1 1], 6, "beta", "0, 2")


            # test error is thrown if parameter input is incorrect
            @test_throws ErrorException continuation(hopf2d, [1 1], 6, "bet", 0:0.01:2)
            @test_throws ErrorException continuation(hopf2d, [1 1], 6, 1, 0:0.01:2)
            @test_throws ErrorException continuation(hopf2d, [1 1], 6, [1], 0:0.01:2)

            ProgressMeter.next!(test)

        end

        @testset verbose = true "Output Tests" begin

            # test if found limit cycle matches the analytical solution
            new_par_values, conditions = continuation(hopf2d, [1 1], 6, "beta", 0:0.01:2)
            @test isapprox(conditions[end, 1:2], [0; 0], atol = 1e-8)
            ProgressMeter.next!(test)


            new_par_values, conditions = continuation(predprey, [1 1], 6, "b", 0.2:0.01:0.3)
            @test isapprox(
                new_par_values[findfirst(isapprox.(conditions, 0.270, atol = 0.001))],
                0.26,
                atol = 0.01,
            )
            ProgressMeter.next!(test)

        end
    end

    @testset verbose = true "Finite Difference" begin

        @testset verbose = true "Input Tests" begin

            # Define allowed inputs
            κ, L, T, mx, mt = 1, 1, 0.5, 10, 1000

            # test error is thrown if kappa inputs are incorrect
            @test_throws ErrorException finite_difference(u_I, 1:0.1:10, L, T, mx, mt)
            @test_throws ErrorException finite_difference(u_I, "1", L, T, mx, mt)
            @test_throws ErrorException finite_difference(u_I, [1], L, T, mx, mt)


            # test error is thrown if L inputs are incorrect
            @test_throws ErrorException finite_difference(u_I, κ, -1, T, mx, mt)
            @test_throws ErrorException finite_difference(u_I, κ, "1", T, mx, mt)
            @test_throws ErrorException finite_difference(u_I, κ, [1], T, mx, mt)


            # test error is thrown if T inputs are incorrect
            @test_throws ErrorException finite_difference(u_I, κ, L, -1, mx, mt)
            @test_throws ErrorException finite_difference(u_I, κ, L, "1", mx, mt)
            @test_throws ErrorException finite_difference(u_I, κ, L, [1], mx, mt)


            # test error is thrown if mx inputs are incorrect
            @test_throws ErrorException finite_difference(u_I, κ, L, T, -10, mt)
            @test_throws ErrorException finite_difference(u_I, κ, L, T, "10", mt)
            @test_throws ErrorException finite_difference(u_I, κ, L, T, [10], mt)


            # test error is thrown if mt inputs are incorrect
            @test_throws ErrorException finite_difference(u_I, κ, L, T, mx, -1000)
            @test_throws ErrorException finite_difference(u_I, κ, L, T, mx, "1000")
            @test_throws ErrorException finite_difference(u_I, κ, L, T, mx, [1000])


            # test error is thrown if method input is incorrect
            @test_throws ErrorException finite_difference(u_I, κ, L, T, mx, mt, method = 1)
            @test_throws ErrorException finite_difference(
                u_I,
                κ,
                L,
                T,
                mx,
                mt,
                method = "pseudo_arclengths",
            )

            ProgressMeter.next!(test)

        end

        @testset verbose = true "Output Tests" begin

            κ, L, T, mx, mt = 1, 1, 0.5, 10, 1000

            # test if predicted behaviour matches the exact solution
            x_est, u_j = finite_difference(u_I, κ, L, T, mx, mt)
            @test isapprox(u_j[1], u_exact(0, T), atol = 1e-10)
            @test isapprox(u_j[end], u_exact(L, T), atol = 1e-10)

            ProgressMeter.next!(test)

            println("\n")

        end
    end
end



"""
Plotting stuff
"""

if save_figures

    ## Plot of ODE for t = 20
    t = 0:0.02:20
    u0 = [-1 0]
    hopf_plot = plot_ode(hopf2d, u0, t, ["u1 (x)" "u2 (y)"])
    savefig(hopf_plot, "./output/Hopf Plot.png")
    savefig(hopf_plot, "./output/Hopf Plot.html")
    ProgressMeter.next!(test)


    ## Plot 3d phase portrait
    hopf_phase_portrait = plot_phase_portrait(hopf2d, u0, t, ["u1 (x)" "u2 (y)"])
    savefig(hopf_phase_portrait, "./output/Hopf Phase Portrait.png")
    savefig(hopf_phase_portrait, "./output/Hopf Phase Portrait.html")
    ProgressMeter.next!(test)


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
    savefig(sol_plot, "./output/Hopf Solution Plot.html")
    ProgressMeter.next!(test)


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
    savefig(fd_method_comp, "./output/Finite Difference Method Comparison.html")
    ProgressMeter.next!(test)
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
    savefig(hopf_phase_portrait_3d, "./output/3d Phase Portrait - Hopf.html")
    ProgressMeter.next!(test)



    ## Plot 3d phase portrait of Cheng-Wang System
    t = 0:0.01:100
    u0 = [1 0 1]
    cheng_wang_phase_portrait_3d = plot_phase_portrait_3d(cheng_wang, u0, t)
    savefig(cheng_wang_phase_portrait_3d, "./output/3d Phase Portrait - Cheng-Wang.png")
    savefig(cheng_wang_phase_portrait_3d, "./output/3d Phase Portrait - Cheng-Wang.html")
    ProgressMeter.next!(test)



    ## Plot 3d phase portrait of Lorenz System
    t = 0:0.01:100
    u0 = [1 1 1]
    lorenz_phase_portrait_3d = plot_phase_portrait_3d(lorenz, u0, t)
    savefig(lorenz_phase_portrait_3d, "./output/3d Phase Portrait - Lorenz.png")
    savefig(lorenz_phase_portrait_3d, "./output/3d Phase Portrait - Lorenz.html")
    ProgressMeter.next!(test)



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
    savefig(fd_method_3d_comp, "./output/Finite Difference Method 3d Comparison.html")

end

println("\n")
ProgressMeter.finish!(test)
