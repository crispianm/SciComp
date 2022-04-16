function euler_step(f, x0, tn, Δt, arg...)
	"""
    Performs a Euler step from t0 to (tn + Δt).

	Parameters:
		f (function): Function which returns a singular value or 1 x n matrix of values.
        x0 (matrix): Matrix of initial values in the 1 x n form, eg: [1] or [1 1].
		tn (float): Starting value for Euler step.
		Δt (float): Step size.
		arg (list, optional): Arguments to pass to f.

	Returns:
        [xn1, tn1] (vector): List of values after step at time tn1.
	"""
    xn1 = x0 + Δt*f(x0, tn, arg...)
    tn1 = tn + Δt

    return [xn1, tn1]
end


function heun3_step(f, x0, tn, Δt, arg...)
	"""
    Performs a Heun 3 step from t0 to (tn + Δt).

	Parameters:
		f (function): Function which returns a singular value or 1 x n matrix of values.
        x0 (matrix): Matrix of initial values in the 1 x n form, eg: [1] or [1 1].
		tn (float): Starting value for Euler step.
		Δt (float): Step size.
		arg (list, optional): Arguments to pass to f.

	Returns:
		[xn1, tn1] (vector): List of values after step at time tn1.
	"""
    k1 = f(x0, tn, arg...)
    k2 = f((x0 + (Δt*k1)/3), (tn + Δt/3), arg...)
    k3 = f((x0 + (2*Δt*k2)/3), (tn + (2*Δt)/3), arg...)

    xn1 = x0 + Δt*(k1 + 3*k3)/4
    tn1 = tn + Δt

    return [xn1, tn1]
end


function ralston4_step(f, x0, tn, Δt, arg...)
	"""
    Performs a Ralston 4 step from t0 to (tn + Δt).

	Parameters:
		f (function): Function which returns a singular value or 1 x n matrix of values.
        x0 (matrix): Matrix of initial values in the 1 x n form, eg: [1] or [1 1].
		tn (float): Starting value for Euler step.
		Δt (float): Step size.
		arg (list, optional): Arguments to pass to f.

	Returns:
		[xn1, tn1] (vector): List of values after step at time tn1.
	"""
    k1 = f(x0, tn, arg...)
    k2 = f((x0 + (Δt*k1*0.4)), (tn + (Δt*0.4)), arg...)
    k3 = f((x0 + Δt*(0.29697761*k1 + 0.15875964*k2)), (tn + Δt*0.45573725), arg...)
    k4 = f((x0 + Δt*(0.21810040*k1 - 3.05096516*k2 + 3.83286476*k3)), (tn + Δt), arg...)
    

    xn1 = x0 + Δt*(0.17476028*k1 - 0.55148066*k2 + 1.20553560*k3 + 0.17118468*k4)
    tn1 = tn + Δt

    return [xn1, tn1]
end


function rk4_step(f, x0, tn, Δt, arg...)
	"""
    Performs a Runge-Kutta 4 step from t0 to (tn + Δt).

	Parameters:
		f (function): Function which returns a singular value or 1 x n matrix of values.
        x0 (matrix): Matrix of initial values in the 1 x n form, eg: [1] or [1 1].
		tn (float): Starting value for Euler step.
		Δt (float): Step size.
		arg (list, optional): Arguments to pass to f.

	Returns:
		[xn1, tn1] (vector): List of values after step at time tn1.
	"""
    k1 = f(x0, tn, arg...)
    k2 = f((x0 + (Δt*k1)/2), (tn + Δt/2), arg...)
    k3 = f((x0 + (Δt*k2)/2), (tn + Δt/2), arg...)
    k4 = f((x0 + Δt*k3), (tn + Δt), arg...)

    xn1 = x0 + Δt*(k1 + 2*k2 + 2*k3 + k4)/6
    tn1 = tn + Δt

    return [xn1, tn1]
end


function three_eight_rule_step(f, x0, tn, Δt, arg...)
	"""
    Performs a 3/8 Rule step from t0 to (tn + Δt).

	Parameters:
		f (function): Function which returns a singular value or 1 x n matrix of values.
        x0 (matrix): Matrix of initial values in the 1 x n form, eg: [1] or [1 1].
		tn (float): Starting value for Euler step.
		Δt (float): Step size.
		arg (list, optional): Arguments to pass to f.

	Returns:
		[xn1, tn1] (vector): List of values after step at time tn1.
	"""
    k1 = f(x0, tn, arg...)
    k2 = f((x0 + (Δt*k1)/3), (tn + Δt/3), arg...)
    k3 = f((x0 + Δt*((-k1/3) + k2)), (tn + 2*Δt/3), arg...)
    k4 = f((x0 + Δt*(k1 - k2 + k3)), (tn + Δt), arg...)

    xn1 = x0 + Δt*(k1 + 3*k2 + 3*k3 + k4)/8
    tn1 = tn + Δt

    return [xn1, tn1]
end


function solve_to(f, x0, t1, t2, Δt, method, arg...)
	"""
    Solves ODE from 

	Parameters:
		f (function): Function which returns a singular value or 1 x n matrix of values.
        x0 (matrix): Matrix of initial values in the 1 x n form, eg: [1] or [1 1].
		t1 (float): Starting time.
        t2 (float): Ending time.
		Δt (float): Step size.
        method (function): method used to approximate solution
            Allowable methods:
                euler_step
                heun3_step
                ralston4_step
                rk4_step
                three_eight_rule_step
		arg (list, optional): Arguments to pass to f.

	Returns:
		x (matrix): 1 x n estimate of ODE at time t2
	"""
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


function solve_ode(f, x0, t, method=rk4_step, Δt=0.001, arg...)
    """
    Solves a provided ODE, f, along time input t with initia condition(s) x0.

	Parameters:
		f (function): Function which returns a singular value or 1 x n matrix of values.
        x0 (matrix): Matrix of initial values in the 1 x n form, eg: [1] or [1 1].
		t (array): Time values to solve between. The first element of t should be the initial value, 0.
        method (function): method used to approximate solution
            Allowable methods:
                euler_step
                heun3_step
                ralston4_step
                rk4_step
                three_eight_rule_step
		Δt (float): Step size.
		arg (list, optional): Arguments to pass to f.

	Returns:
        x_series: List of solutions for x at each time value in t
	"""

    # Error handling

    if t[1] != 0
        throw(error("Please make sure the first value of the time series is 0."))
    elseif !isa(x0, Array)
        throw(error("Please make sure the initial condition is a 1 x n matrix."))
    end

    allowable_inputs = [euler_step, heun3_step, ralston4_step, rk4_step, three_eight_rule_step]

    if method ∉ allowable_inputs
        println("Method not assigned, please enter either:
        euler_step, heun3_step, ralston4_step, rk4_step, or three_eight_rule_step\n")
        println("Defaulting to rk4_step.")
        method = rk4_step
    else
        # println("Method successfully set as ", method)
    end

    # Computation

    x_series = Matrix{Float64}(undef, 0, length(x0))
    x_series = [x_series; x0]

    x = x0
    for i = 1:(length(t)-1)
        x = solve_to(f, x, t[i], t[i + 1], Δt, method)
        x_series = [x_series; x]
    end
    
    return x_series
end


