function euler_step(f, x0, tn, h, arg...)

    xn1 = x0 + h*f(x0, tn, arg...)
    tn1 = tn + h

    return [xn1, tn1]
end


function rk4_step(f, x0, tn, h, arg...)

    k1 = f(x0, tn, arg...)
    k2 = f((x0 .+ (h.*k1)./2), (tn .+ h./2), arg...)
    k3 = f((x0 .+ (h.*k2)./2), (tn .+ h./2), arg...)
    k4 = f((x0 .+ h.*k3), (tn .+ h), arg...)

    xn1 = x0 + h*(k1 .+ 2*k2 .+ 2*k3 .+ k4)./6
    tn1 = tn + h

    return [xn1, tn1]
end


function solve_to(f, x0, t1, t2, deltat_max, method, arg...)

    timesteps = floor((t2 - t1) / deltat_max)

    x = x0
    t = t1

    for i = 1:timesteps
        x, t = method(f, x, t, deltat_max, arg...)
    end

    if t != t2
        h = t2 - t
        x, t = method(f, x, t, h, arg...)
    end

    return x
end


function solve_ode(f, x0, t, method, deltat_max, arg...)
    
    if t[1] != 0
        throw(error("Please make sure the first value of the time series is 0."))
    end
    
    x_series = Matrix{Float64}(undef, 0, length(x0))
    x_series = [x_series; x0]

    x = x0
    for i = 1:(length(t)-1)
        x = solve_to(f, x, t[i], t[i + 1], deltat_max, method)
        x_series = [x_series; x]
    end
    
    return x_series
end


