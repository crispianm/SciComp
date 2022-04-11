include("../visualisation.jl")
include("../ode_solver.jl")
include("../numerical_shooting.jl")

function hopf3d(u, beta, sigma=-1.0)

    u1, u2, u3 = u
    du1dt = beta*u1 - u2 + sigma*u1*(u1^2 + u2^2)
    du2dt = u1 + beta*u2 + sigma*u2*(u1^2 + u2^2)
    du3dt = -u3

    return [du1dt du2dt du3dt]
end


function hopf_sol3d(t, beta, theta=0.0)

    u1 = √(beta) * cos.(t .+ theta)
    u2 = √(beta) * sin.(t .+ theta)
    u3 = exp.(-t)

    return u1, u2, u3
end

t = 0:0.1:100
u0 = [1 1 1]

plot_phase_portrait_3d(hopf3d, [1 0 -2], t, ["u1" "u2" "u3"])