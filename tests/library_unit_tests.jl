# Import Subsystems
include("../examples/example_functions.jl")
include("../finite_difference.jl")
include("../numerical_continuation.jl")
include("../numerical_shooting.jl")
include("../ode_solver.jl")
using Test
using ProgressMeter


"""
    Test script to run the ODE solver for a variety of differential equations
    and check them against explicit solutions.

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


println("\n")

ProgressMeter.finish!(test)
