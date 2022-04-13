function forward_euler(init_condition, kappa, L, T, mx, mt)
    
    # Set up the numerical environment variables
    x = LinRange(0, L, mx+1)     # mesh points in space
    t = LinRange(0, T, mt+1)     # mesh points in time
    deltax = x[2] - x[1]            # gridspacing in x
    deltat = t[2] - t[1]            # gridspacing in t
    lmbda = kappa*deltat/(deltax^2)    # mesh fourier number
    println("deltax = ", deltax)
    println("deltat = ", deltat)
    println("lambda = ", lmbda)

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
        for i in 2:mx
            u_jp1[i] = u_j[i] + lmbda*(u_j[i-1] - 2*u_j[i] + u_j[i+1])
        end
        # Boundary conditions
        u_jp1[1] = 0; u_jp1[mx+1] = 0
            
        # Save u_j at time t[j+1]
        u_j[:] = u_jp1[:]
    end


    return x, u_j
end