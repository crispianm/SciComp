function G(f, u0, t0, T, arg...)

    include("ode_solver.jl")

    F = solve_ode(f, u0, [t0 T], rk4_step, 0.001, arg...)
    diff = u0 .- F[[end],:]

    return diff
end

function shoot(U, f, phase_index=1, arg...)

    u0 = U[:, 1:end-1]
    T = U[end]
    g_estimate = G(f, u0, 0, T, arg...)
    phase_conditions = f(u0, 0, arg...)[phase_index]
   
    return [g_estimate phase_conditions]
end