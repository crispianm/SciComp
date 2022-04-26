using NLsolve
using LinearAlgebra
include("./NumericalShooting.jl")


function np_continuation(f, u0, T, par_values, discretisation; arg...)

    """
    Continuation method using the last found solution as an initial guess.
    First solution is found using the inputted initial conditions, u0.

        Parameters:
            f (function): Function which returns a singular value or 1 x n matrix of values.
            u0 (matrix): Matrix of initial values in the 1 x n form, eg: [1] or [1 1].
            T (float): Initial guess for the period.
            par_values (range): Parameter values to solve between, made with colons, eg: 0:0.1:2.
            discretisation (anon. function): The discretisation to use, defined in the "continuation" function below.
            arg (list, optional): Arguments to pass to f.

        Returns:
            par_values, conditions: the parameter values and corresponding solutions.
    """

    conditions = nlsolve((u) -> discretisation(u, par_values[1]), u0).zero
    for parameter in par_values[2:end]
        # Solve using previous initial condition
        x = nlsolve((u) -> discretisation(u, parameter), conditions[[end], :]).zero
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

function pseudo_arclength(f, u0, T, par_values, discretisation; arg...)

    """
    Continuation method using the pseudo arclength equation.

        Parameters:
            f (function): Function which returns a singular value or 1 x n matrix of values.
            u0 (matrix): Matrix of initial values in the 1 x n form, eg: [1] or [1 1].
            T (float): Initial guess for the period.
            par_values (range): Parameter values to solve between, made with colons, eg: 0:0.1:2.
            discretisation (anon. function): The discretisation to use, defined in the "continuation" function below.
            arg (list, optional): Arguments to pass to f.

        Returns:
            par_values, conditions: the parameter values and corresponding solutions.
    """

    # Create a vector of new parameter values to try
    new_par_values = [par_values[1]; par_values[2]]

    # Anonymous function to check if the parameter difference is in the range of par_values
    is_in_range = (value) -> xor(value < par_values[end], value < par_values[1])

    # Use the first solution as an initial guess for the next
    conditions = nlsolve((u) -> discretisation(u, par_values[1]), u0).zero
    conditions =
        [conditions; nlsolve((u) -> discretisation(u, new_par_values[2]), conditions).zero]

    i = 1
    while is_in_range(new_par_values[end])

        # Define augmented state vectors
        v0 = [new_par_values[i] conditions[[i], :]]
        v1 = [new_par_values[i+1] conditions[[i + 1], :]]

        # Find the secant
        secant = v1 - v0

        # Find the pseudo-arclength estimate
        v_pred = v1 + secant
        sol =
            nlsolve(
                (v2) -> [discretisation(v2[:, 2:end], v2[1]) pseudo_arclength_eq(
                    secant,
                    v2,
                    v_pred,
                )],
                v_pred,
            ).zero

        # Append the new condition and parameter values
        conditions = [conditions; sol[:, 2:end]]
        push!(new_par_values, sol[1])

        i += 1
    end

    return new_par_values, conditions
end


function continuation(
    f,
    u0,
    T,
    parameter,
    par_values;
    method = "pseudo_arclength",
    discretisation = "shooting",
    arg...,
)

    """
    Finds a function's solution for each a range of parameter values.

        Parameters:
            f (function): Function which returns a singular value or 1 x n matrix of values.
            u0 (matrix): Matrix of initial values in the 1 x n form, eg: [1] or [1 1].
            T (float): Initial guess for the period.
            parameter (string): The parameter in the system to vary.
                Allowable parameter inputs:
                    "a", "alpha", "b", "beta", "c", "d", "delta", "sigma"
            par_values (range): Parameter values to solve between, made with colons, eg: 0:0.1:2.
            method (string): Method to use to find the solution.
                Allowable method inputs:
                    "pseudo_arclength", "pa", "natural_parameter_continuation", or "npc"
            discretisation (string): The discretisation to use.
                Allowable discretisation inputs:
                    "shooting" or "none"
            arg (list, optional): Arguments to pass to f.

        Returns:
            par_values, conditions: the parameter values and corresponding solutions.

        Example Usage:
            continuation(hopf2d, [1 1], 6, "beta", 0:0.01:2, method="pseudo_arclength", discretisation="shooting")
    """

    ## Error handling
    if !isa(u0, Array)
        error(
            "Please make sure the initial condition is a 1 x n matrix.\neg: [1] or [1 1].",
        )
    elseif size(u0)[1] != 1
        error(
            "Please make sure the initial condition is a 1 x n matrix.\neg: [1] or [1 1].",
        )
    elseif !isa(parameter, String)
        error("Please enter a string for the phase index.")
    elseif !isa(discretisation, String)
        error("Please enter a string for the discretisation.")
    elseif !isa(method, String)
        error("Please enter a string for the method.")
    elseif typeof(T) ∉ (Float64, Int)
        error("Please enter a single integer or float for the period.")
    end

    # Check par_values is a range.
    try
        step(par_values)
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
    allowable_discretisations = ["shooting", "none"]
    if string(discretisation) ∉ allowable_discretisations
        error("Discretisation not assigned, please enter either:
        shooting or none.")
    end

    # Check method is allowed
    allowable_methods = ["pseudo_arclength", "pa", "natural_parameter_continuation", "npc"]
    if string(method) ∉ allowable_methods
        error("Method not assigned, please enter either:
        pseudo_arclength or natural_parameter_continuation.")
    end

    # convert any potential ints to floats for use in nlsolve
    u0 = [u0 T]
    u0 = [Float64(number) for number in u0]
    par_values = [Float64(value) for value in par_values]

    # Handle parameter and discretisation assignment
    if discretisation == "shooting"

        if parameter == "a"
            discretisation = (u0, par) -> shoot(f, u0, phase_index = 0, a = par)
        elseif parameter == "alpha"
            discretisation = (u0, par) -> shoot(f, u0, phase_index = 0, alpha = par)
        elseif parameter == "b"
            discretisation = (u0, par) -> shoot(f, u0, phase_index = 0, b = par)
        elseif parameter == "beta"
            discretisation = (u0, par) -> shoot(f, u0, phase_index = 0, beta = par)
        elseif parameter == "c"
            discretisation = (u0, par) -> shoot(f, u0, phase_index = 0, c = par)
        elseif parameter == "d"
            discretisation = (u0, par) -> shoot(f, u0, phase_index = 0, d = par)
        elseif parameter == "delta"
            discretisation = (u0, par) -> shoot(f, u0, phase_index = 0, delta = par)
        elseif parameter == "sigma"
            discretisation = (u0, par) -> shoot(f, u0, phase_index = 0, sigma = par)
        end

    elseif discretisation == "none"

        if parameter == "a"
            discretisation = (u0, par) -> f(u0, a = par)
        elseif parameter == "alpha"
            discretisation = (u0, par) -> f(u0, alpha = par)
        elseif parameter == "b"
            discretisation = (u0, par) -> f(u0, b = par)
        elseif parameter == "beta"
            discretisation = (u0, par) -> f(u0, beta = par)
        elseif parameter == "c"
            discretisation = (u0, par) -> f(u0, c = par)
        elseif parameter == "d"
            discretisation = (u0, par) -> f(u0, d = par)
        elseif parameter == "delta"
            discretisation = (u0, par) -> f(u0, delta = par)
        elseif parameter == "sigma"
            discretisation = (u0, par) -> f(u0, sigma = par)
        end

    else
        error("Invalid discretisation.")
    end


    ## Computation

    if method == "natural_parameter_continuation" || method == "npc"
        new_par_values, conditions = np_continuation(f, u0, T, par_values, discretisation)
    elseif method == "pseudo_arclength" || method == "pa"
        new_par_values, conditions = pseudo_arclength(f, u0, T, par_values, discretisation)
    else
        error("Unknown method: ", method)
    end

    return new_par_values, conditions
end
