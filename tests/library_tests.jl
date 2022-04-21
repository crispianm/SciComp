using PlotlyJS
using Test
include("../ode_solver.jl")
include("../numerical_shooting.jl")
include("../examples/example_functions.jl")
include("../visualisation.jl")


"""
Test script to run the ODE solver for a variety of differential equations and check them against explicit solutions.


Set the save_figures and save_figures_3d variables to true to save pngs in the output folder.
"""
save_figures = true
save_figures_3d = true


"""
Run the tests!
"""



# x0 = [1]
# t = [0, 1]

# estimate = solve_ode(f, x0, t, "rk4", 0.001)[end][1]

# @testset "Systems Tests" begin

#     @testset "Hopf Week 17" begin
#         @test foo("cat") == 9
#         @test foo("dog") == foo("cat")
#     end

#     @testset "Arrays $i" for i in 1:3
#         @test foo(zeros(i)) == i^2
#         @test foo(ones(i)) == i^2
#     end

# end









"""
Plotting stuff
"""

if save_figures

    ## Plot of ODE for t = 20
    t = 0:0.02:20
    u0 = [-1 0]
    hopf_plot = plot_ode(hopf2d, u0, t, ["u1 (x)" "u2 (y)"])
    savefig(hopf_plot, "./output/Hopf Plot.png")


    ## Plot phase portrait
    hopf_phase_portrait = plot_phase_portrait(hopf2d, u0, t, ["u1 (x)" "u2 (y)"])
    savefig(hopf_phase_portrait, "./output/Hopf Phase Portrait.png")


    ## Plot ODE analytic and numerical solutions to see if they match
    # Find limit cycle
    u0, T = find_limit_cycle(hopf2d, [-1 0], 6.28)
    println("u0: ", u0)
    println("Period: ", T)

    t = 0:0.01:T
    hopf_solution = hopf2d_sol(t, beta=1, theta=pi)

    # Numerical solution
    numerical_sol = solve_ode(hopf2d, u0, t, "rk4")

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
end



"""
3d plotting stuff

This is just to check if the visualisation and the numerical shooting emthods are working,
    by visually comparing output to that of known solutions. And I found it kinda fun.
"""

if save_figures_3d

    # Hopf 3d System
    t = 0:0.1:200
    u0 = [1 0 -2]
    hopf_phase_portrait_3d = plot_phase_portrait_3d(hopf3d, u0, t);
    savefig(hopf_phase_portrait_3d, "./output/3d Phase Portrait - Hopf.png")

    # Cheng-Wang System
    t = 0:0.01:100
    u0 = [1 0 1]
    cheng_wang_phase_portrait_3d = plot_phase_portrait_3d(cheng_wang, u0, t);
    savefig(cheng_wang_phase_portrait_3d, "./output/3d Phase Portrait - Cheng-Wang.png")

    # Lorenz System
    t = 0:0.01:100
    u0 = [1 1 1]
    lorenz_phase_portrait_3d = plot_phase_portrait_3d(lorenz, u0, t);
    savefig(lorenz_phase_portrait_3d, "./output/3d Phase Portrait - Lorenz.png")

end