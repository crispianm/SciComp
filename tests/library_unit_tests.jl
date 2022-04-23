using PlotlyJS
using Test
include("../ode_solver.jl")
include("../numerical_continuation.jl")
include("../numerical_shooting.jl")
include("../examples/example_functions.jl")
include("../visualisation.jl")
include("../finite_difference.jl")


#  Test script to run the ODE solver for a variety of differential equations
#  and check them against explicit solutions.
 
#  Set the save_figures and save_figures_3d variables to true
#  to save pngs in the output folder.


save_figures = false
save_figures_3d = false


"""
Run the tests!
"""

test = ProgressUnknown("Performing tests: ", spinner=true)

@testset verbose = true "System Unit Tests" begin

    @testset verbose = true "ode_solver" begin

        @testset verbose = true "Input Tests" begin

            # test error is thrown if t=0 is not included in t
            t = 0:0.1:1
            bad_t = 1:0.1:2
            @test_throws ErrorException solve_ode(f2, [1], bad_t, method="rk4")


            # test error is thrown if x0 is not a matrix
            @test_throws ErrorException solve_ode(f2, 1, t, method="rk4")
            @test_throws ErrorException solve_ode(f2, [1;0], t, method="rk4")
            @test_throws ErrorException solve_ode(f2, [1, 0], t, method="rk4")

        
            # test error is thrown if x0 is not a 1xn matrix
            @test_throws ErrorException solve_ode(f2, [1 2; 3 4], t, method="rk4")


            # test error is thrown if x0 is not correct length
            @test_throws ErrorException solve_ode(f2, [1], t, method="rk4")

            ProgressMeter.next!(test)

        end

        @testset verbose = true "Output Tests" begin

            # test if solve_ode estimates a simple ODE correctly
            x0 = [1]
            t = 0:1
            e_estimate = solve_ode(f, x0, t, method="rk4")[end][1]
            @test isapprox(e_estimate, ℯ)


            # test if solve_ode estimates systems of ODEs correctly
            x0 = [1 0]
            t = 0:0.1:1
            f2_sol = f2_solution(x0, t)
            f2_numerical_sol = solve_ode(f2, x0, t, method="rk4")
            @test  all(isapprox.(f2_numerical_sol, f2_sol, atol=1e-6))

            ProgressMeter.next!(test)

        end
    end

    @testset verbose = true "numerical_shooting" begin

        @testset verbose = true "Input Tests" begin

            # test error is thrown if T is not an integer or float
            T = 10
            bad_T = 1:0.1:2
            @test_throws ErrorException find_limit_cycle(f2, [1], bad_T)
        
        
            # test error is thrown if u0 is not a matrix
            @test_throws ErrorException find_limit_cycle(f2, 1, T)
            @test_throws ErrorException find_limit_cycle(f2, [1;0], T)
            @test_throws ErrorException find_limit_cycle(f2, [1, 0], T)
        
        
            # test error is thrown if u0 is not a 1xn matrix
            @test_throws ErrorException find_limit_cycle(f2, [1 2; 3 4], T)
        
        
            # test error is thrown if u0 is not correct length
            @test_throws ErrorException find_limit_cycle(f2, [1], 10)
        
        
            # test error is thrown if phase_index is a positive integer
            @test_throws ErrorException find_limit_cycle(f2, [1 1], 10, phase_index=-1)
            @test_throws ErrorException find_limit_cycle(f2, [1 1], 10, phase_index=[3])
            @test_throws ErrorException find_limit_cycle(f2, [1 1], 10, phase_index=2.5)
            @test_throws ErrorException find_limit_cycle(f2, [1 1], 10, phase_index="0")

            ProgressMeter.next!(test)
        
        end

        @testset verbose = true "Output Tests" begin

            # test if found limit cycle matches the analytical solution
            u0, T = find_limit_cycle(hopf2d, [-1 0], 6)
            @test isapprox(T, 2*pi)


            # test if solve_ode estimates a Hopf ODE correctly
            t = 0:0.1:T
            hopf_solution = hopf2d_sol(t, beta=1, theta=pi) # adjusted for phase
            hopf_numerical_sol = solve_ode(hopf2d, u0, t, method="rk4")
            @test  all(isapprox.(hopf_numerical_sol, hopf_solution, atol=1e-6))


            # test if arguments are being passed to the function correctly
            u0_0, T = find_limit_cycle(hopf2d, [-1 0], 6)
            u0_1, T = find_limit_cycle(hopf2d, [-1 0], 6, beta = 2)
            u0_3, T = find_limit_cycle(hopf2d, [-1 0], 6, beta = 2, sigma=-1.2)
            u0_2, T = find_limit_cycle(hopf2d, [-1 0], 6, sigma=-1.2)
            @test u0_0 != u0_1
            @test u0_0 != u0_2
            @test u0_0 != u0_3

            ProgressMeter.next!(test)

        end

        @testset verbose = true "Higher Dimensional Output" begin

            # test if found limit cycle matches the analytical solution
            u0, T = find_limit_cycle(hopf3d, [1 1 1], 6)
            @test isapprox(T, 2*pi)


            u0, T = find_limit_cycle(hopf3d, [1 1 1], 10)
            hopf_solution = hopf3d_sol(u0, 0:0.1:T; theta=0) # adjusted for phase
            hopf_numerical_sol = solve_ode(hopf3d, u0, 0:0.1:T, method="rk4")
            @test all(isapprox.(hopf_numerical_sol, hopf_solution, atol=1e-5))

            ProgressMeter.next!(test)
            
        end
    end

    @testset verbose = true "numerical_continuation" begin

        @testset verbose = true "Input Tests" begin
    
            # test error is thrown if T is not an integer or float
            T = 10
            bad_T = 1:0.1:2
            @test_throws ErrorException continuation(hopf2d, [1 1], bad_T, "beta", 0:0.01:2)
        
        
            # test error is thrown if u0 is not a matrix
            @test_throws ErrorException continuation(hopf2d, 1, 6, "beta", 0:0.01:2)
            @test_throws ErrorException continuation(hopf2d, [1;0], 6, "beta", 0:0.01:2)
            @test_throws ErrorException continuation(hopf2d, [1, 0], 6, "beta", 0:0.01:2)
        
        
            # test error is thrown if u0 is not a 1xn matrix
            @test_throws ErrorException continuation(hopf2d, [1 2; 3 4], 6, "beta", 0:0.01:2)
        
        
            # test error is thrown if u0 is not correct length
            @test_throws ErrorException continuation(hopf2d, [1], 6, "beta", 0:0.01:2)
        
            
            # test error is thrown if method input is incorrect
            @test_throws ErrorException continuation(hopf2d, [1 1], 6, "beta", 0:0.01:2,
                method = 1)
            @test_throws ErrorException continuation(hopf2d, [1 1], 6, "beta", 0:0.01:2,
                method = "pseudo_arclengths")
    
            
            # test error is thrown if discretisation input is incorrect
            @test_throws ErrorException continuation(hopf2d, [1 1], 6, "beta", 0:0.01:2,
                discretisation = 2)
            @test_throws ErrorException continuation(hopf2d, [1 1], 6, "beta", 0:0.01:2,
                discretisation = "shoot")
    
    
            # test error is thrown if par_values input is incorrect
            @test_throws ErrorException continuation(hopf2d, [1 1], 6, "beta", [0 2])
            @test_throws ErrorException continuation(hopf2d, [1 1], 6, "beta", [0.1 0.5 1.0 1.5 2.0])
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
            @test isapprox(conditions[199, 1:2], [0; 0], atol=1e-8)


            new_par_values, conditions = continuation(predprey, [1 1], 6, "b", 0.2:0.01:0.3)
            @test isapprox(new_par_values[findfirst(isapprox.(conditions, 0.270, atol=0.001))], 0.26, atol=0.01)
            
            ProgressMeter.next!(test)
    
        end
    end

    @testset verbose = true "finite_difference" begin

        @testset verbose = true "Input Tests" begin

            # test error is thrown if t=0 is not included in t
            t = 0:0.1:1
            bad_t = 1:0.1:2
            @test_throws ErrorException solve_ode(f2, [1], bad_t, method="rk4")

            ProgressMeter.next!(test)

        end

        @testset verbose = true "Output Tests" begin
            
            # test if found limit cycle matches the analytical solution
            u0, T = find_limit_cycle(hopf2d, [-1 0], 6)
            @test isapprox(T, 2*pi)

            ProgressMeter.next!(test)

            println("\n")

        end
    end
end

println("\n")
ProgressMeter.finish!(test)












"""
Plotting stuff
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
    hopf_solution = hopf2d_sol(t, beta=1, theta=pi)

    # Numerical solution
    numerical_sol = solve_ode(hopf2d, u0, t)

    x_plot = scatter(
        x=t, y=numerical_sol[:,1],
        mode="lines", name="u1 (x) - numerical",
        showlegend=true)
    y_plot = scatter(
        x=t, y=numerical_sol[:,2],
        mode="lines", name="u2 (y) - numerical",
        showlegend=true)

    # Analytic solution
    u1_plot = scatter(x=t, y=hopf_solution[:,1], mode="lines", name="u1 (x)  - analytic", showlegend=true)
    u2_plot = scatter(x=t, y=hopf_solution[:,2], mode="lines", name="u2 (y)  - analytic", showlegend=true)

    layout = Layout(xaxis_title = "time")

    sol_plot = plot([x_plot, y_plot, u1_plot, u2_plot], layout)
    savefig(sol_plot, "./output/Hopf Solution Plot.png")


    ## Plot finite_difference comparison plot
    
end



"""
3d plotting stuff

This is just to check if the visualisation and the numerical shooting emthods are working,
    by visually comparing output to that of known solutions. And I found it kinda fun.
"""

if save_figures_3d

    ## Plot 3d phase portrait of Hopf 3d System
    t = 0:0.1:200
    u0 = [1 0 -2]
    hopf_phase_portrait_3d = plot_phase_portrait_3d(hopf3d, u0, t);
    savefig(hopf_phase_portrait_3d, "./output/3d Phase Portrait - Hopf.png")

    ## Plot 3d phase portrait of Cheng-Wang System
    t = 0:0.01:100
    u0 = [1 0 1]
    cheng_wang_phase_portrait_3d = plot_phase_portrait_3d(cheng_wang, u0, t);
    savefig(cheng_wang_phase_portrait_3d, "./output/3d Phase Portrait - Cheng-Wang.png")

    ## Plot 3d phase portrait of Lorenz System
    t = 0:0.01:100
    u0 = [1 1 1]
    lorenz_phase_portrait_3d = plot_phase_portrait_3d(lorenz, u0, t);
    savefig(lorenz_phase_portrait_3d, "./output/3d Phase Portrait - Lorenz.png")

end