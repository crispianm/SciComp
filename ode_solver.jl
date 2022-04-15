function euler_step(f, x0, tn, Δt, arg...)

    xn1 = x0 + Δt*f(x0, tn, arg...)
    tn1 = tn + Δt

    return [xn1; tn1]
end


function rk4_step(f, x0, tn, Δt, arg...)

    k1 = f(x0, tn, arg...)
    k2 = f((x0 + (Δt*k1)/2), (tn + Δt/2), arg...)
    k3 = f((x0 + (Δt*k2)/2), (tn + Δt/2), arg...)
    k4 = f((x0 + Δt*k3), (tn + Δt), arg...)

    xn1 = x0 + Δt*(k1 + 2*k2 + 2*k3 + k4)/6
    tn1 = tn + Δt

    return [xn1; tn1]
end


function solve_to(f, x0, t1, t2, Δt, method, arg...)

    timesteps = floor((t2 - t1) / Δt)

    x = x0
    t = t1

    for i = 1:timesteps
        x, t = method(f, x, t, Δt, arg...)
    end

    if t != t2
        Δt = t2 - t
        x, t = method(f, x, t, Δt, arg...)
    end

    return x
end


function solve_ode(f, x0, t, method, Δt, arg...)
    
    if t[1] != 0
        throw(error("Please make sure the first value of the time series is 0."))
    end
    
    x_series = Matrix{Float64}(undef, 0, length(x0))
    x_series = [x_series; x0]

    x = x0
    for i = 1:(length(t)-1)
        x = solve_to(f, x, t[i], t[i + 1], Δt, method)
        x_series = [x_series; x]
    end
    
    return x_series
end


