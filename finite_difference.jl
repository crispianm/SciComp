using PlotlyJS
using LinearAlgebra

function forward_euler(u_I, lmbda, L, T, mx, mt, x, arg...)

    # Create A_FE matrix
    A_FE = Tridiagonal(ones(mx-1)*(lmbda), ones(mx)*(1 - 2*lmbda), ones(mx-1)*(lmbda))

    # Set up the solution variables
    u_j = zeros(size(x))        # u at current time step
    u_jp1 = zeros(size(x))      # u at next time step

    # Set initial condition
    for i in 1:mx+1
        u_j[i] = u_I(x[i])
    end

    # Solve the PDE: loop over all time points
    for j in 1:mt
        # Forward Euler timestep at inner mesh points
        # PDE discretised at position x[i], time t[j]

        u_jp1[2:end] = A_FE * u_j[2:end]

        # Boundary conditions
        u_jp1[1] = 0
        u_jp1[mx+1] = 0
            
        # Save u_j at time t[j+1]
        u_j[:] = u_jp1[:]
    end

    return x, u_j
end

function backward_euler(u_I, lmbda, L, T, mx, mt, x, arg...)

    # Create A_BE matrix
    A_BE = Tridiagonal(ones(mx-1)*(-lmbda), ones(mx)*(1 + 2*lmbda), ones(mx-1)*(-lmbda))
    println("X = ", mx)
    # Set up the solution variables
    u_j = zeros(size(x))        # u at current time step
    u_jp1 = zeros(size(x))      # u at next time step

    # Set initial condition
    for i in 1:mx+1
        u_j[i] = u_I(x[i])
    end

    # Solve the PDE: loop over all time points
    for j in 1:mt
        # Forward Euler timestep at inner mesh points
        # PDE discretised at position x[i], time t[j]

        u_jp1[2:end] = A_BE \ u_j[2:end]

        # Boundary conditions
        u_jp1[1] = 0
        u_jp1[mx+1] = 0
            
        # Save u_j at time t[j+1]
        u_j[:] = u_jp1[:]
    end

    return x, u_j
end

function crank_nicholson(u_I, lmbda, L, T, mx, mt, x, arg...)

    # Create A_CN and B_CN matrices
    A_CN = Tridiagonal(ones(mx-1)*(-lmbda/2), ones(mx)*(1 + lmbda), ones(mx-1)*(-lmbda/2))
    B_CN = Tridiagonal(ones(mx-1)*(lmbda/2), ones(mx)*(1 - lmbda), ones(mx-1)*(lmbda/2))

    # Set up the solution variables
    u_j = zeros(size(x))        # u at current time step
    u_jp1 = zeros(size(x))      # u at next time step

    # Set initial condition
    for i in 1:mx+1
        u_j[i] = u_I(x[i])
    end

    # Solve the PDE: loop over all time points
    for j in 1:mt
        # Forward Euler timestep at inner mesh points
        # PDE discretised at position x[i], time t[j]

        u_jp1[2:end] = A_CN \ B_CN * u_j[2:end]

        # Boundary conditions
        u_jp1[1] = 0
        u_jp1[mx+1] = 0
            
        # Save u_j at time t[j+1]
        u_j[:] = u_jp1[:]
    end

    return x, u_j
end


function finite_difference(u_I, kappa, L, T, mx, mt, method, arg...)

    # Set up the numerical environment variables
    x = 0:(L/mx):L     # mesh points in space
    t = 0:(T/mt):T     # mesh points in time
    deltax = x[2] - x[1]            # gridspacing in x
    deltat = t[2] - t[1]            # gridspacing in t
    lmbda = kappa*deltat/(deltax^2)    # mesh fourier number

    if method == "forward_euler" || method == "fe"
        x1, u_j = forward_euler(u_I, lmbda, L, T, mx, mt, x)
    elseif method == "backward_euler" || method == "be"
        x1, u_j = backward_euler(u_I, lmbda, L, T, mx, mt, x)
    elseif method == "crank_nicholson" || method == "cn"
        x1, u_j = crank_nicholson(u_I, lmbda, L, T, mx, mt, x)
    else
        error("Unknown method: ", method)
    end

    return x1, u_j
end