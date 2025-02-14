# Various ODE functions and their solutions are included in this file.

"""
dx/dt = x
"""

function f(u, t; arg...)

    if size(u)[end] != 1
        error("Please enter only one initial condition for the function.")
    end

    return u
end

function f_solution(u, t; arg...)

    if size(u)[end] != 1
        error("Please enter only one initial condition for the function.")
    end

    return u * exp.(t)
end



"""
System of equations definitions for d^2x/dt^2 = -x, equivalent to the system of equations
    dx/dt = y
and
    dy/dt = -x.
"""

function f2(u, t; arg...)

    if size(u)[end] != 2
        error("Please enter two initial conditions for the function.")
    end

    u1, u2 = u

    du1dt = u2
    du2dt = -u1
    
    return [du1dt du2dt]
end

function f2_solution(u, t; arg...)

    c1, c2 = u

    x = c1 * cos.(t) .+ c2 * sin.(t)
    y = c2 * cos.(t) .- c1 * sin.(t)

    return [x y]
end



"""
System of equations definitions for d^2x/dt^2 = -x, equivalent to the system of equations
    dx/dt = y
and
    dy/dt = -x.
"""

function predprey(u, t; a = 1, b = 0.2, d = 0.1, arg...)

    if size(u)[end] != 2
        error("Please enter two initial conditions for the function.")
    end

    x, y = u

    x_dot = x .* (1 .- x) .- (a .* x .* y) ./ (d .+ x)
    y_dot = b .* y .* (1 .- (y ./ x))

    return [x_dot y_dot]
end



"""
Algebraic equation for numerical continuation exercises
"""

function algebraic(x; c = 1, arg...)

    return x .^ 3 .- x .+ c
end



"""
Hopf bifurcation for numerical continuation exercises
"""

function hopf2d(u, t; beta = 1, sigma = -1.0, arg...)

    if size(u)[end] != 2
        error("Please enter two initial conditions for the function.")
    end

    u1, u2 = u
    du1dt = beta .* u1 .- u2 .+ sigma .* u1 .* (u1 .^ 2 + u2 .^ 2)
    du2dt = u1 .+ beta .* u2 .+ sigma .* u2 .* (u1 .^ 2 + u2 .^ 2)

    return [du1dt du2dt]
end

function hopf2d_modified(u, t; beta = 1, sigma = -1.0, arg...)

    if size(u)[end] != 2
        error("Please enter two initial conditions for the function.")
    end

    u1, u2 = u
    du1dt =
        beta .* u1 .- u2 .- sigma .* u1 .* (u1 .^ 2 .+ u2 .^ 2) .+
        sigma .* u1 .* (u1 .^ 2 .+ u2 .^ 2) .^ 2
    du2dt =
        u1 .+ beta .* u2 .- sigma .* u2 .* (u1 .^ 2 .+ u2 .^ 2) .+
        sigma .* u2 .* (u1 .^ 2 .+ u2 .^ 2) .^ 2

    return [du1dt du2dt]
end

function hopf2d_sol(t; beta = 1, theta = 0.0, arg...)

    u1 = √(beta) .* cos.(t .+ theta)
    u2 = √(beta) .* sin.(t .+ theta)

    return [u1 u2]
end



"""
Hopf bifurcation in 3d
"""

function hopf3d(u, t; beta = 1, sigma = -1.0, arg...)

    if size(u)[end] != 3
        error("Please enter three initial conditions for the function.")
    end

    u1, u2, u3 = u
    du1dt = beta .* u1 .- u2 .+ sigma .* u1 .* (u1 .^ 2 .+ u2 .^ 2)
    du2dt = u1 .+ beta .* u2 .+ sigma .* u2 .* (u1 .^ 2 .+ u2 .^ 2)
    du3dt = -u3

    return [du1dt du2dt du3dt]
end

function hopf3d_sol(u, t; beta = 1, theta = 0.0, arg...)

    c1, c2, c3 = u

    u1 = √(beta) .* cos.(t .+ theta)
    u2 = √(beta) .* sin.(t .+ theta)
    u3 = c3 * exp.(-t) # correct?

    return [u1 u2 u3]
end



"""
Finite Difference Functions
"""

# Example initial temperature distribution
function u_I(x; L = 1.0, arg...)

    y = @. sin.(pi * x / L)

    return y
end

# The exact solution of the above function
function u_exact(x, t; L = 1.0, κ = 1.0, arg...)

    y = @. exp(-κ * ((pi / L)^2) * t) * sin(pi * x / L)

    return y
end



"""
Other Functions (fun to plot)
"""

function cheng_wang(u, t; a = -0.01, arg...)

    """
    Good initial condition: [1 0 1]
    """

    x, y, z = u
    dxdt = y
    dydt = z
    dzdt = a .- y .- x .^ 2 .- x .* z .+ 3 * y .^ 2

    return [dxdt dydt dzdt]
end

function lorenz(u, t; beta = (8 / 3), sigma = 10.0, rho = 28.0, arg...)

    """
    Good initial condition: [1 1 1]
    """

    x, y, z = u
    dxdt = sigma .* (y .- x)
    dydt = x * (rho - z) - y
    dzdt = x * y - beta .* z

    return [dxdt dydt dzdt]
end

function van_der_pol(u, t; μ = 1, arg...)

    """
    Good initial condition: [2 0]
    """

    x, y = u
    dxdt = y
    dydt = μ .* (1 .- x .^ 2) .* y .- x

    return [dxdt dydt]
end
