include("../visualisation.jl")
include("../ode_solver.jl")
include("../numerical_shooting.jl")

function hopf3d(u, t, beta=1, sigma=-1.0, arg...)

    u1, u2, u3 = u
    du1dt = beta*u1 - u2 + sigma*u1*(u1^2 + u2^2)
    du2dt = u1 + beta*u2 + sigma*u2*(u1^2 + u2^2)
    du3dt = -u3

    return [du1dt du2dt du3dt]
end

function hopf_sol3d(t, beta=1, theta=0.0, arg...)

    u1 = √(beta) * cos.(t .+ theta)
    u2 = √(beta) * sin.(t .+ theta)
    u3 = exp.(-t)

    return u1, u2, u3
end

function cheng_wang(u, t, a=-0.01, arg...)

    x, y, z = u
    dxdt = y
    dydt = z
    dzdt = a - y - x^2 - x*z + 3*y^2

    return [dxdt dydt dzdt]
end

function lorenz(u, t, beta=(8/3), sigma=10.0, rho=28.0, arg...)

    x, y, z = u
    dxdt = sigma*(y - x)
    dydt = x*(rho - z) - y
    dzdt = x*y - beta*z

    return [dxdt dydt dzdt]
end


t = 0:0.01:100

plot_phase_portrait_3d(hopf3d, [1 0 -2], t)
plot_phase_portrait_3d(cheng_wang, [1 0 1], t)
plot_phase_portrait_3d(lorenz, [1 1 1], t)
