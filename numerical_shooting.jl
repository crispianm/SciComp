using NLsolve
include("ode_solver.jl")

function G(f, u0, t0, T, arg...)

    F = solve_ode(f, u0, [t0 T], rk4_step, 0.001, arg...)
    g = u0 .- F[[end], :]
    return g
end

function shoot(f, u, phase_index=1, arg...)

    u0 = u[:, 1:end .!= end]
    T = u[end]
    
    g_est = G(f, u0, 0, T, arg...)
    f_at_index = f(u0, 0, arg...)[phase_index]
    
    return [g_est f_at_index]
end


function find_limit_cycle(f, u0, T, phase_index=1, arg...)

    U = [u0 T]
    solution = nlsolve((u) -> shoot(f, u), U).zero

    u0_1 = solution[:, 1:end-1]
    T_1 = solution[end]

    return u0_1, T_1
end