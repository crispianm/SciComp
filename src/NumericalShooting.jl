using NLsolve
include("ODESolver.jl")

function G(f, u0, T; arg...)

    """
    Solves the ODE, f, and returns u0 minus the solution.

        Parameters:
            f (function): Function which returns a singular value or 1 x n matrix of values.
            u0 (matrix): Matrix of initial values in the 1 x n form, eg: [1] or [1 1].
            T (float): Period estimate
            arg (list, optional): Arguments to pass to f.

        Returns:
            u0 minus the solution of ODE f.
    """

    F = solve_ode(f, u0, [0 T]; arg...)
    g = u0 .- F[[end], :]

    return g
end

function shoot(f, u; phase_index = 0, arg...)

    """
    Defines the system to be solved using in find_limit_cycle.

        Parameters:
            f (function): Function which returns a singular value or 1 x n matrix of values.
            u (matrix): Matrix of initial value(s) and time period in the 1 x n form, eg: [1 20] or [1 1 20].
            phase_index (int, optional): Variable for the phase condition, eg: d/dt(phase_index) = 0
            arg (list, optional): Arguments to pass to f.

        Returns:
            A matrix of G's estimate and the phase condition.

        Example Usage:
            shoot(hopf2d, [1 1 20], phase_index = 0)
    """

    u0 = u[:, 1:end.!=end]
    T = u[end]

    if T == 0
        error("T cannot be 0.")
    end

    G_estimate = G(f, u0, T; arg...)
    phase_condition = f(u0, 0; arg...)[phase_index+1]

    return [G_estimate phase_condition]
end


function find_limit_cycle(f, u0, T; phase_index = 0, arg...)

    """
    Finds u0 and T such that the phase condition of f is satisfied.
        
        Parameters:
            f (function): Function which returns a singular value or 1 x n matrix of values.
            u0 (matrix): Initial guess for the initial conditions in the 1 x n form, eg: [1] or [1 1].
            T (float): Initial guess for the period.
            phase_index (int, optional): Variable for the phase condition, eg: d/dt(phase_index) = 0
            arg (list, optional): Arguments to pass to f.
        
        Returns:
            found_u0, found_T, the initial conditions and period of the limit cycle.

        Example usage:
            find_limit_cycle(hopf2d, [-1 0], 20, phase_index = 0)
    """

    # Error handling
    if !isa(u0, Array)
        error(
            "Please make sure the initial condition is a 1 x n matrix.\neg: [1] or [1 1].",
        )
    elseif size(u0)[1] != 1
        error(string(
            "Please make sure the initial condition is a size 1 x n matrix, instead of a size ", size(u0), " ",typeof(u0),
            "\neg: [1] or [1 1]."
        ))
    elseif !isa(phase_index, Int)
        error("Please enter an integer for the phase index.")
    elseif phase_index < 0
        error("Please enter a positive integer for the phase index.")
    elseif typeof(T) ∉ (Float64, Int)
        error("Please enter a single integer or float for the period.")
    end

    U = [u0 T]
    U = [Float64(number) for number in U] # convert ints to floats for use in nlsolve

    solution = nlsolve((u) -> shoot(f, u; phase_index, arg...), U).zero

    found_u0 = solution[:, 1:end-1]
    found_T = solution[end]

    return found_u0, found_T
end
