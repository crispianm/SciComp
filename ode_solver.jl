# function euler_step(f, x0, tn, h, arg...)

#     euler = Any[]
#     xn1 = x0 + h*f(x0, tn, arg...)
#     # print(xn1)
#     tn1 = tn + h

#     push!(euler, xn1)
#     push!(euler, tn1)

#     return euler
# end

function euler_step(f, x0, tn, h, arg...)

    xn1 = x0 + h*f(x0, tn, arg...)
    tn1 = tn + h

    return [xn1; tn1]
end

function rk4_step(f, x0, tn, h, arg...)

    rk4 = Any[]

    k1 = f(x0, tn, arg...)
    k2 = f((x0 .+ (h.*k1)./2), (tn .+ h./2), arg...)
    k3 = f((x0 .+ (h.*k2)./2), (tn .+ h./2), arg...)
    k4 = f((x0 .+ h.*k3), (tn .+ h), arg...)

    xn1 = x0 + h*(k1 .+ 2*k2 .+ 2*k3 .+ k4)./6
    tn1 = tn + h

    push!(rk4, xn1)
    push!(rk4, tn1)

    return rk4
end


function solve_to(f, x0, t1, t2, method, deltat_max, arg...)

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
    
    if(!isapprox(t[1], 0.0; atol=eps(Float64), rtol=0))
        throw(error("Please make sure the first value of the time series is 0."))
    end
    
    x_series = Any[]
    push!(x_series, x0[1])

    x = x0
    for i = 1:(length(t)-1)
        x = solve_to(f, x, t[i], t[i + 1], method, deltat_max, arg...)
        push!(x_series, x)
    end

    return x_series
end


