function G(f, u0, t0=0, T, arg...)

    include("ode_solver.jl")

    F = solve_ode(f, u0, [t0 T], rk4_step, 0.001, arg...)
    g = u0 .- F[[end],:]

    return g
end

function shoot(f, u0, T, phase_index=1, arg...)

    g_est = G(f, u0, 0, T, arg...)
    f_at_index = f(u0, 0, arg...)[phase_index]
   
    return [g_est f_at_index]
end