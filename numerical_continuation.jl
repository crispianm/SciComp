function np_continuation(f, x0, T, parameter, par_values; discretisation, arg...)
    
    """
    Attempts to find a function's solution for each parameter value using the last found solution as an initial guess.
    First solution is found using the inputted initial conditions, x0.

        Parameters:
            f (function): Function which returns a singular value or 1 x n matrix of values.
                The parameter 
            x0 (matrix): Matrix of initial values in the 1 x n form, eg: [1] or [1 1].
            T (float): Initial guess for the period.
            parameter (string): The parameter in the system to vary.
                Allowable method inputs: a, alpha, b, beta, c, d, delta, sigma
            par_values (range): Parameter values to solve between, made with colons, eg: 0:0.1:2.
            discretisation (string): The discretisation to use, either "shooting" or "none."
            arg (list, optional): Arguments to pass to f.

        Returns:
            par_values, conditions: the parameter values and corresponding solutions.

        Example Usage:
            np_continuation(hopf2d, [1 1], 6, "beta", 0:0.01:2, discretisation="shooting")
    """

    ## Error handling
    if !isa(x0, Array)
        error("Please make sure the initial condition is a 1 x n matrix.\neg: [1] or [1 1].")
    elseif size(x0)[1] != 1
        error("Please make sure the initial condition is a 1 x n matrix.\neg: [1] or [1 1].")
    elseif !isa(parameter, String)
        error("Please enter a string for the phase index.")
    elseif !isa(discretisation, String)
        error("Please enter a string for the discretisation.")
    elseif typeof(T) ∉ (Float64, Int)
        error("Please enter a single integer or float for the period.")
    end

    # Check par_values is a range.
    try 
        step(par_values);
    catch 
        error("Please enter a range for the parameter values.")
    end

    # Check parameter is allowed
    allowable_parameters = ["a", "b", "c", "d", "alpha", "beta", "delta", "sigma"]
    if string(parameter) ∉ allowable_parameters
        error("Parameter not assigned, please enter either:
        a, b, c, d, alpha, beta, delta, or sigma")
    end

    # Check discretisation method is allowed
    allowable_parameters = ["shooting", "none"]
    if string(discretisation) ∉ allowable_parameters
        error("Discretisation not assigned, please enter either:
        shooting or none.")
    end

    # convert ints to floats for use in nlsolve
    x0 = [x0 T]
    x0 = [Float64(number) for number in x0] 
    
    # Handle parameter and discretisation assignment
    if discretisation == "shooting"

        if parameter == "a"
            discretisation = (u0, par) -> shoot(f, u0, phase_index=0, a=par)
        elseif parameter == "alpha"
            discretisation = (u0, par) -> shoot(f, u0, phase_index=0, alpha=par)
        elseif parameter == "b"
            discretisation = (u0, par) -> shoot(f, u0, phase_index=0, b=par)
        elseif parameter == "beta"
            discretisation = (u0, par) -> shoot(f, u0, phase_index=0, beta=par)
        elseif parameter == "c"
            discretisation = (u0, par) -> shoot(f, u0, phase_index=0, c=par)
        elseif parameter == "d"
            discretisation = (u0, par) -> shoot(f, u0, phase_index=0, d=par)
        elseif parameter == "delta"
            discretisation = (u0, par) -> shoot(f, u0, phase_index=0, delta=par)
        elseif parameter == "sigma"
            discretisation = (u0, par) -> shoot(f, u0, phase_index=0, sigma=par)
        end

    elseif discretisation == "none"
        
        if parameter == "a"
            discretisation = (u0, par) -> f(u0, a=par)
        elseif parameter == "alpha"
            discretisation = (u0, par) -> f(u0, alpha=par)
        elseif parameter == "b"
            discretisation = (u0, par) -> f(u0, b=par)
        elseif parameter == "beta"
            discretisation = (u0, par) -> f(u0, beta=par)
        elseif parameter == "c"
            discretisation = (u0, par) -> f(u0, c=par)
        elseif parameter == "d"
            discretisation = (u0, par) -> f(u0, d=par)
        elseif parameter == "delta"
            discretisation = (u0, par) -> f(u0, delta=par)
        elseif parameter == "sigma"
            discretisation = (u0, par) -> f(u0, sigma=par)
        end
       
    else
        error("Invalid discretisation.")
    end


    ## Computation
    conditions = nlsolve((u) -> discretisation(u, par_values[1]), x0).zero
    for parameter in par_values[2:end]
        # solve using previous initial condition
        x = nlsolve((u) -> discretisation(u, parameter), conditions[[end],:]).zero
        conditions = [conditions; x]
    end

    return par_values, conditions
end


function pseudo_arclength_eq(secant, v, v_pred)

    """
    Finds the pseudo arclength estimate.

        Parameters:
            secant (matrix): Matrix of secant values.
            v (matrix): Matrix of values at the secant points.
            v_pred (matrix): Matrix of values at the secant points, predicted by a continuation method.

        Returns:
            The pseudo arclength estimate.

    """

    estimate = dot(secant, v - v_pred)

    return estimate
end


function pseudo_arclength(f, x0, T, parameter, par_values; discretisation="shooting", arg...)

    """
    Attempts to find a function's solution using the pseudo arclength method

        Parameters:
            f (function): Function which returns a singular value or 1 x n matrix of values.
                The parameter 
            x0 (matrix): Matrix of initial values in the 1 x n form, eg: [1] or [1 1].
            T (float): Initial guess for the period.
            parameter (string): The parameter in the system to vary.
                Allowable method inputs: a, alpha, b, beta, c, d, delta, sigma
            par_values (range): Parameter values to solve between, made with colons, eg: 0:0.1:2.
            discretisation (string): The discretisation to use, either "shooting" or "none."
            arg (list, optional): Arguments to pass to f.

        Returns:
            new_par_values, conditions: the parameter values and corresponding solutions.

        Example Usage:
            pseudo_arclength(hopf2d, [1 1], 6, "beta", 0:0.01:2, discretisation="shooting")
    """

    ## Error handling
    if !isa(x0, Array)
        error("Please make sure the initial condition is a 1 x n matrix.\neg: [1] or [1 1].")
    elseif size(x0)[1] != 1
        error("Please make sure the initial condition is a 1 x n matrix.\neg: [1] or [1 1].")
    elseif !isa(parameter, String)
        error("Please enter a string for the phase index.")
    elseif !isa(discretisation, String)
        error("Please enter a string for the discretisation.")
    elseif typeof(T) ∉ (Float64, Int)
        error("Please enter a single integer or float for the period.")
    end

    # Check par_values is a range.
    try 
        step(par_values);
    catch 
        error("Please enter a range for the parameter values.")
    end

    # Check parameter is allowed
    allowable_parameters = ["a", "b", "c", "d", "alpha", "beta", "delta", "sigma"]
    if string(parameter) ∉ allowable_parameters
        error("Parameter not assigned, please enter either:
        a, b, c, d, alpha, beta, delta, or sigma")
    end

    # Check discretisation method is allowed
    allowable_parameters = ["shooting", "none"]
    if string(discretisation) ∉ allowable_parameters
        error("Discretisation not assigned, please enter either:
        shooting or none.")
    end

    # convert ints to floats for use in nlsolve
    x0 = [x0 T]
    x0 = [Float64(number) for number in x0] 
    
    # Handle parameter and discretisation assignment
    if discretisation == "shooting"
        if parameter == "a"
            discretisation = (u0, par) -> shoot(f, u0, phase_index=0, a=par)
        elseif parameter == "alpha"
            discretisation = (u0, par) -> shoot(f, u0, phase_index=0, alpha=par)
        elseif parameter == "b"
            discretisation = (u0, par) -> shoot(f, u0, phase_index=0, b=par)
        elseif parameter == "beta"
            discretisation = (u0, par) -> shoot(f, u0, phase_index=0, beta=par)
        elseif parameter == "c"
            discretisation = (u0, par) -> shoot(f, u0, phase_index=0, c=par)
        elseif parameter == "d"
            discretisation = (u0, par) -> shoot(f, u0, phase_index=0, d=par)
        elseif parameter == "delta"
            discretisation = (u0, par) -> shoot(f, u0, phase_index=0, delta=par)
        elseif parameter == "sigma"
            discretisation = (u0, par) -> shoot(f, u0, phase_index=0, sigma=par)
        end
    elseif discretisation == "none" 
        if parameter == "a"
            discretisation = (u0, par) -> f(u0, a=par)
        elseif parameter == "alpha"
            discretisation = (u0, par) -> f(u0, alpha=par)
        elseif parameter == "b"
            discretisation = (u0, par) -> f(u0, b=par)
        elseif parameter == "beta"
            discretisation = (u0, par) -> f(u0, beta=par)
        elseif parameter == "c"
            discretisation = (u0, par) -> f(u0, c=par)
        elseif parameter == "d"
            discretisation = (u0, par) -> f(u0, d=par)
        elseif parameter == "delta"
            discretisation = (u0, par) -> f(u0, delta=par)
        elseif parameter == "sigma"
            discretisation = (u0, par) -> f(u0, sigma=par)
        end
    else
        error("Invalid discretisation.")
    end


    ## Computation

    # Create a list of new parameter values to try
    new_par_values = [par_values[1]; par_values[2]]

    # Check if the parameter difference is positive or negative 
    if par_values[end] - par_values[1] < 0
        end_function = (value) -> value > par_values[end]
    else
        end_function = (value) -> value < par_values[end]
    end

    # Use the first solution as an initial guess for the next
    conditions = nlsolve((u) -> discretisation(u, par_values[1]), x0).zero
    conditions = [conditions; nlsolve((u) -> discretisation(u, new_par_values[2]), conditions).zero]
   
    i = 1
    while end_function(new_par_values[end])

        # Define augmented state vectors
        v0 = [new_par_values[i] conditions[[i],:]]
        v1 = [new_par_values[i+1] conditions[[i+1],:]]
        
        # Find the secant
        secant = v1 - v0
        
        # Find the pseudo-arclength estimate
        v_pred = v1 + secant
        sol = nlsolve((v2) -> [discretisation(v2[:,2:end], v2[1]) pseudo_arclength_eq(secant, v2, v_pred)], v_pred).zero
        
        # Append the new condition and parameter value
        conditions = [conditions; sol[:,2:end]]
        push!(new_par_values, sol[1])

        i += 1
    end
    
    return new_par_values, conditions
end