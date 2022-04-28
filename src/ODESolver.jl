function euler_step(ode, u0, tn, Δt; arg...)

    """
       Performs a Euler step from t0 to (tn + Δt).

    Parameters:
    	ode (function): Function which returns a singular value or 1 x n matrix of values.
        u0 (matrix): Matrix of initial values in the 1 x n form, eg: [1] or [1 1].
    	tn (float): Starting value for Euler step.
    	Δt (float): Step size.
    	arg (list, optional): Arguments to pass to f.

    Returns:
           [xn1, tn1] (vector): List of values after step at time tn1.
    """

    xn1 = u0 + Δt * ode(u0, tn; arg...)
    tn1 = tn + Δt

    return [xn1, tn1]
end

function heun3_step(ode, u0, tn, Δt; arg...)

    """
       Performs a Heun 3 step from t0 to (tn + Δt).

    Parameters:
    	ode (function): Function which returns a singular value or 1 x n matrix of values.
        u0 (matrix): Matrix of initial values in the 1 x n form, eg: [1] or [1 1].
    	tn (float): Starting value for Euler step.
    	Δt (float): Step size.
    	arg (list, optional): Arguments to pass to f.

    Returns:
    	[xn1, tn1] (vector): List of values after step at time tn1.
    """

    k1 = ode(u0, tn; arg...)
    k2 = ode((u0 + (Δt * k1) / 3), (tn + Δt / 3); arg...)
    k3 = ode((u0 + (2 * Δt * k2) / 3), (tn + (2 * Δt) / 3); arg...)

    xn1 = u0 + Δt * (k1 + 3 * k3) / 4
    tn1 = tn + Δt

    return [xn1, tn1]
end

function ralston4_step(ode, u0, tn, Δt; arg...)

    """
       Performs a Ralston 4 step from t0 to (tn + Δt).

    Parameters:
    	ode (function): Function which returns a singular value or 1 x n matrix of values.
        u0 (matrix): Matrix of initial values in the 1 x n form, eg: [1] or [1 1].
    	tn (float): Starting value for Euler step.
    	Δt (float): Step size.
    	arg (list, optional): Arguments to pass to f.

    Returns:
    	[xn1, tn1] (vector): List of values after step at time tn1.
    """

    k1 = ode(u0, tn; arg...)
    k2 = ode((u0 + (Δt * k1 * 0.4)), (tn + (Δt * 0.4)); arg...)
    k3 = ode((u0 + Δt * (0.29697761 * k1 + 0.15875964 * k2)), (tn + Δt * 0.45573725); arg...)
    k4 = ode(
        (u0 + Δt * (0.21810040 * k1 - 3.05096516 * k2 + 3.83286476 * k3)),
        (tn + Δt);
        arg...,
    )


    xn1 = u0 + Δt * (0.17476028 * k1 - 0.55148066 * k2 + 1.20553560 * k3 + 0.17118468 * k4)
    tn1 = tn + Δt

    return [xn1, tn1]
end

function rk4_step(ode, u0, tn, Δt; arg...)

    """
       Performs a Runge-Kutta 4 step from t0 to (tn + Δt).

    Parameters:
    	ode (function): Function which returns a singular value or 1 x n matrix of values.
        u0 (matrix): Matrix of initial values in the 1 x n form, eg: [1] or [1 1].
    	tn (float): Starting value for Euler step.
    	Δt (float): Step size.
    	arg (list, optional): Arguments to pass to f.

    Returns:
    	[xn1, tn1] (vector): List of values after step at time tn1.
    """

    k1 = ode(u0, tn; arg...)
    k2 = ode((u0 + (Δt * k1) / 2), (tn + Δt / 2); arg...)
    k3 = ode((u0 + (Δt * k2) / 2), (tn + Δt / 2); arg...)
    k4 = ode((u0 + Δt * k3), (tn + Δt); arg...)

    xn1 = u0 + Δt * (k1 + 2 * k2 + 2 * k3 + k4) / 6
    tn1 = tn + Δt

    return [xn1, tn1]
end

function three_eighths_rule_step(ode, u0, tn, Δt; arg...)

    """
       Performs a 3/8 Rule step from t0 to (tn + Δt).

    Parameters:
    	ode (function): Function which returns a singular value or 1 x n matrix of values.
        u0 (matrix): Matrix of initial values in the 1 x n form, eg: [1] or [1 1].
    	tn (float): Starting value for Euler step.
    	Δt (float): Step size.
    	arg (list, optional): Arguments to pass to f.

    Returns:
    	[xn1, tn1] (vector): List of values after step at time tn1.
    """

    k1 = ode(u0, tn; arg...)
    k2 = ode((u0 + (Δt * k1) / 3), (tn + Δt / 3); arg...)
    k3 = ode((u0 + Δt * ((-k1 / 3) + k2)), (tn + 2 * Δt / 3); arg...)
    k4 = ode((u0 + Δt * (k1 - k2 + k3)), (tn + Δt); arg...)

    xn1 = u0 + Δt * (k1 + 3 * k2 + 3 * k3 + k4) / 8
    tn1 = tn + Δt

    return [xn1, tn1]
end

function solve_to(ode, u0, t1, t2, Δt, method; arg...)

    """
       Solves ODE from 

    Parameters:
    	ode (function): Function which returns a singular value or 1 x n matrix of values.
        u0 (matrix): Matrix of initial values in the 1 x n form, eg: [1] or [1 1].
    	t1 (float): Starting time.
        t2 (float): Ending time.
    	Δt (float): Step size.
           method (function): method used to approximate solution
               Allowable methods:
                   euler_step
                   heun3_step
                   ralston4_step
                   rk4_step
                   three_eighths_rule_step
    	arg (list, optional): Arguments to pass to f.

    Returns:
    	x (matrix): 1 x n estimate of ODE at time t2
    """

    timesteps = floor((t2 - t1) / Δt)

    x = u0
    t = t1

    for i = 1:timesteps
        x, t = method(ode, x, t, Δt; arg...)
    end

    if t != t2
        Δt = t2 - t
        x, t = method(ode, x, t, Δt; arg...)
    end

    return x
end

function solve_ode(ode, u0, t; method = "rk4", Δt = 0.01, arg...)

    """
    Solves a provided ODE, ode, along time input t with initial condition(s) u0.
    
    Parameters:
        ode (function): Function which returns a singular value or 1 x n matrix of values.
        u0 (matrix): Matrix of initial values in the 1 x n form, eg: [1] or [1 1].
        t (array or range): Time values to solve between. The first element of t should be the initial value, 0.
        method (string): method used to approximate solution
            Allowable method inputs:
                "euler", "heun3", "ralston4", "rk4", "three_eighths_rule"
        Δt (float): Step size.
        arg (list, optional): Arguments to pass to f.
    
    Returns:
        x_series: List of solutions for x at each time value in t
    """

    # Error handling
    if t[1] != 0
        error("Please make sure the first value of the time series is 0.")
    elseif length(t) < 2
        error(string("Please make sure the time series has at least 2 values, instead of ", length(t),"."))
    elseif !isa(u0, Array)
        error(string(
            "Please make sure the initial condition is a 1 x n matrix, not a ", typeof(u0), 
            ".\neg: [1] or [1 1].",
        ))
    elseif size(u0)[1] != 1
        error(string(
            "Please make sure the initial condition is a 1 x n matrix, instead of ", size(u0),
            "\neg: [1] or [1 1]."
        ))
    end

    # Check method is allowed
    allowable_inputs = ["euler", "heun3", "ralston4", "rk4", "three_eighths_rule"]

    if string(method) ∉ allowable_inputs
        println("Method not assigned, please enter either:
        euler, heun3, ralston4, rk4, or three_eighths_rule")
        println("Defaulting to rk4_step.")
        method = rk4_step
    else
        method = string(method, "_step")
        method = getfield(Main, Symbol(method))
        # println("Method successfully set as ", method)
    end

    # Computation
    x_series = Matrix{Float64}(undef, 0, length(u0))
    x_series = [x_series; u0]

    x = u0
    for i = 1:(length(t)-1)
        x = solve_to(ode, x, t[i], t[i+1], Δt, method; arg...)
        x_series = [x_series; x]
    end

    return x_series
end
