using LinearAlgebra

function zero_boundary(x, t, arg...)

    """
    Boundary condition = 0.

        Parameters:
            x: x-coordinate
            t: time

        Returns:
            0: Always. It's very useful.
    """

    return 0
end

function forward_euler(u_initial, λ, mx, mt, x, t; boundary = zero_boundary, arg...)

    """
    Forward Euler estimate of PDE, for use in the finite_difference function

        Parameters:
            u_initial (function): Initial temperature distribution function.
            λ (float): Lambda parameter, found inside the finite_difference function.
            mx (int): Number of gridpoints in space.
            mt (int): Number of gridpoints in time.
            x (range): Range of x values, found inside the finite_difference function.
            t (range): Range of time values, found inside the finite_difference function.
            boundary (function): Boundary condition function. Defaults to zero_boundary
            arg (list, optional): Arguments to pass to u_initial.

        Returns:
            x_est, u_j: The values of u at each x over time T.
    """

    # Create A_FE matrix
    A_FE = Tridiagonal(ones(mx - 1) * (λ), ones(mx) * (1 - 2 * λ), ones(mx - 1) * (λ))

    # Set up the solution variables
    u_j = zeros(size(x))        # u at current time step
    u_jp1 = zeros(size(x))      # u at next time step

    # Set initial condition
    for i = 1:mx+1
        u_j[i] = u_initial(x[i], L = x[end], arg...)
    end

    # Solve the PDE: loop over all time points
    for j = 1:mt
        # Forward Euler timestep at inner mesh points
        # PDE discretised at position x[i], time t[j]
        u_jp1[2:end] = A_FE * u_j[2:end]

        # Boundary conditions
        u_jp1[1] = boundary(0, t[j], arg...)
        u_jp1[mx+1] = boundary(x[end], t[j], arg...)

        # Save u_j at time t[j+1]
        u_j = u_jp1
    end

    return x, u_j
end

function backward_euler(u_initial, λ, mx, mt, x, t; boundary = zero_boundary, arg...)

    """
    Backward Euler estimate of PDE, for use in the finite_difference function

        Parameters:
            u_initial (function): Initial temperature distribution function.
            λ (float): Lambda parameter, found inside the finite_difference function.
            mx (int): Number of gridpoints in space.
            mt (int): Number of gridpoints in time.
            x (range): Range of x values, found inside the finite_difference function.
            t (range): Range of time values, found inside the finite_difference function.
            boundary (function): Boundary condition function. Defaults to zero_boundary
            arg (list, optional): Arguments to pass to u_initial.

        Returns:
            x_est, u_j: The values of u at each x over time T.
    """

    # Create A_BE matrix
    A_BE = Tridiagonal(ones(mx - 1) * (-λ), ones(mx) * (1 + 2 * λ), ones(mx - 1) * (-λ))

    # Set up the solution variables
    u_j = zeros(size(x))        # u at current time step
    u_jp1 = zeros(size(x))      # u at next time step
    # Set initial condition
    for i = 1:mx+1
        u_j[i] = u_initial(x[i], L = x[end], arg...)
    end

    # Solve the PDE: loop over all time points
    for j = 1:mt
        # Forward Euler timestep at inner mesh points
        # PDE discretised at position x[i], time t[j]
        u_jp1[2:end] = A_BE \ u_j[2:end]

        # Boundary conditions
        u_jp1[1] = boundary(0, t[j], arg...)
        u_jp1[mx+1] = boundary(x[end], t[j], arg...)

        # Save u_j at time t[j+1]
        u_j = u_jp1
    end

    return x, u_j
end

function crank_nicholson(u_initial, λ, mx, mt, x, t; boundary = zero_boundary, arg...)

    """
    Crank Nicholson estimate of PDE, for use in the finite_difference function

        Parameters:
            u_initial (function): Initial temperature distribution function.
            λ (float): Lambda parameter, found inside the finite_difference function.
            mx (int): Number of gridpoints in space.
            mt (int): Number of gridpoints in time.
            x (range): Range of x values, found inside the finite_difference function.
            t (range): Range of time values, found inside the finite_difference function.
            boundary (function): Boundary condition function. Defaults to zero_boundary
            arg (list, optional): Arguments to pass to u_initial.

        Returns:
            x_est, u_j: The values of u at each x over time T.
    """

    # Create A_CN and B_CN matrices
    A_CN = Tridiagonal(ones(mx - 1) * (-λ / 2), ones(mx) * (1 + λ), ones(mx - 1) * (-λ / 2))
    B_CN = Tridiagonal(ones(mx - 1) * (λ / 2), ones(mx) * (1 - λ), ones(mx - 1) * (λ / 2))

    # Set up the solution variables
    u_j = zeros(size(x))        # u at current time step
    u_jp1 = zeros(size(x))      # u at next time step

    # Set initial condition
    for i = 1:mx+1
        u_j[i] = u_initial(x[i], L = x[end], arg...)
    end

    # Solve the PDE: loop over all time points
    for j = 1:mt
        # Forward Euler timestep at inner mesh points
        # PDE discretised at position x[i], time t[j]
        u_jp1[2:end] = A_CN \ B_CN * u_j[2:end]

        # Boundary conditions
        u_jp1[1] = boundary(0, t[j], arg...)
        u_jp1[mx+1] = boundary(x[end], t[j], arg...)

        # Save u_j at time t[j+1]
        u_j = u_jp1
    end

    return x, u_j
end


function finite_difference(
    u_initial,
    κ,
    L,
    T,
    mx,
    mt;
    boundary = zero_boundary,
    method = "cn",
    arg...,
)

    """
    Solves PDE using finite differences.

        Parameters:
            u_initial (function): Initial temperature distribution function.
            κ (float): Diffusion constant.
            L (float): Length of spatial domain.
            T (float): Total time to solve for.
            mx (int): Number of gridpoints in space.
            mt (int): Number of gridpoints in time.
            boundary (function): Boundary condition function. Defaults to zero_boundary.
            method (str): The finite difference method to use. Defaults to "cn"
                Allowable method inputs:
                    "forward_euler", "fe", "backward_euler", "be", "crank_nicholson", or "cn"
            arg (list, optional): Arguments to pass to u_initial.

        Returns:
            x_est, u_j: The values of u at each x over time T.
    """

    ## Error handling
    if typeof(κ) ∉ (Int, Float64)
        error("Please enter a single number for κ.")
    elseif typeof(L) ∉ (Int, Float64)
        error("Please enter a single number for L.")
    elseif L < 0
        error("Please enter a positive number for L.")
    elseif typeof(T) ∉ (Int, Float64)
        error("Please enter a single number for T.")
    elseif T < 0
        error("Please enter a positive number for T.")
    elseif typeof(mx) != Int
        error("Please enter a single integer for mx.")
    elseif mx < 0
        error("Please enter a positive integer for mx.")
    elseif typeof(mt) != Int
        error("Please enter a single integer for mt.")
    elseif mt < 0
        error("Please enter a positive integer for mt.")
    elseif typeof(method) != String
        error("Please enter a string for the method. Allowable methods:
        forward_euler, fe, backward_euler, be, crank_nicholson, or cn")
    end

    allowable_methods =
        ["forward_euler", "fe", "backward_euler", "be", "crank_nicholson", "cn"]
    if string(method) ∉ allowable_methods
        error("Method not assigned, please enter either:
        forward_euler, fe, backward_euler, be, crank_nicholson, or cn")
    end

    ## Computation

    # Set up the numerical environment variables
    x = 0:(L/mx):L     # mesh points in space
    t = 0:(T/mt):T     # mesh points in time
    Δx = L / mx          # gridspacing in x
    Δt = T / mt          # gridspacing in t
    λ = κ * Δt / (Δx^2)    # mesh fourier number

    # Assign method and find x_est and u_j
    if method == "forward_euler" || method == "fe"
        x_est, u_j = forward_euler(u_initial, λ, mx, mt, x, t, boundary = boundary, arg...)
    elseif method == "backward_euler" || method == "be"
        x_est, u_j = backward_euler(u_initial, λ, mx, mt, x, t, boundary = boundary, arg...)
    elseif method == "crank_nicholson" || method == "cn"
        x_est, u_j =
            crank_nicholson(u_initial, λ, mx, mt, x, t, boundary = boundary, arg...)
    else
        error("Unknown method: ", method)
    end

    return x_est, u_j
end
