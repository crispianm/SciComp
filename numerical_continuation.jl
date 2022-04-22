function np_continuation(f, x0, T, parameter, par_values; discretisation, arg...)
    
    """
    Attempts to find a function's solution for each parameter value using the last found solution as an initial guess.
    First solution is found using the inputted initial conditions, x0.

        Parameters:
            f (function): Function which returns a singular value or 1 x n matrix of values.
                The parameter 
            x0 (matrix): Matrix of initial values in the 1 x n form, eg: [1] or [1 1].
            T (float): Initial guess for the period.
            parameter: The parameter in the system to vary.
                Allowable method inputs: a, alpha, b, beta, c, d, delta, sigma
            par_values (array or range): Parameter values to solve between.
            discretisation (string): The discretisation to use, either "shooting" or "none."
            arg (list, optional): Arguments to pass to f.

        Returns:
            par_values, conditions: the parameter values and corresponding solutions.

        Example Usage:
            np_continuation(hopf2d, [1 1], 6, 0:0.01:2, discretisation="shooting")
    """

    x0 = [x0 T]
    x0 = [Float64(number) for number in x0] # convert ints to floats for use in nlsolve
    

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


    # Computation
    conditions = nlsolve((u) -> discretisation(u, par_values[1]), x0).zero
    for parameter in par_values[2:end]
        # solve using previous initial condition
        x = nlsolve((u) -> discretisation(u, parameter), conditions[[end],:]).zero
        conditions = [conditions; x]
    end

    return par_values, conditions
end

function pseudo_arclength(f, x0, T, par_values; discretisation, arg...)

    

end
