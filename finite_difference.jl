using PlotlyJS
using LinearAlgebra

function forward_euler(u_I, kappa, L, T, mx, mt)

    # Set up the numerical environment variables
    x = LinRange(0, L, mx+1)     # mesh points in space
    t = LinRange(0, T, mt+1)     # mesh points in time
    deltax = x[2] - x[1]            # gridspacing in x
    deltat = t[2] - t[1]            # gridspacing in t
    lmbda = kappa*deltat/(deltax^2)    # mesh fourier number
    println("deltax = ", deltax)
    println("deltat = ", deltat)
    println("lambda = ", lmbda)

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

    # println("x: ", x)
    # println("u_j: ", u_j) 
    return x, u_j
end

function backward_euler(u_I, kappa, L, T, mx, mt)

    # Set up the numerical environment variables
    x = LinRange(0, L, mx+1)     # mesh points in space
    t = LinRange(0, T, mt+1)     # mesh points in time
    deltax = x[2] - x[1]            # gridspacing in x
    deltat = t[2] - t[1]            # gridspacing in t
    lmbda = kappa*deltat/(deltax^2)    # mesh fourier number
    println("deltax = ", deltax)
    println("deltat = ", deltat)
    println("lambda = ", lmbda)

    # Create A_FE matrix
    A_FE = Tridiagonal(ones(mx-1)*(-lmbda), ones(mx)*(1 + 2*lmbda), ones(mx-1)*(-lmbda))

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

        u_jp1[2:end] = A_FE \ u_j[2:end]

        # Boundary conditions
        u_jp1[1] = 0
        u_jp1[mx+1] = 0
            
        # Save u_j at time t[j+1]
        u_j[:] = u_jp1[:]
    end

    # println("x: ", x)
    # println("u_j: ", u_j) 
    return x, u_j
end

kappa = 1   # diffusion constant
L = 1.0         # length of spatial domain
T = 0.5         # total time to solve for

function u_I(x)
    # initial temperature distribution
    y = sin.(pi*x/L)
    return y
end

function u_exact(x,t)
    # the exact solution
    y = exp.(-kappa*(pi^2/L^2)*t)*sin.(pi*x/L)
    return y
end

# Set numerical parameters
mx = 10     # number of gridpoints in space
mt = 1000   # number of gridpoints in time

# forward euler Estimate
x, u_j = forward_euler(u_I, kappa, L, T, mx, mt);
# Create trace
f_euler = scatter(x=x, y=u_j, mode="markers", name="forward euler", showlegend=true)

# forward euler Estimate
x, u_j = backward_euler(u_I, kappa, L, T, mx, mt);
# Create trace
b_euler = scatter(x=x, y=u_j, mode="markers", name="backward euler", showlegend=true)




# Plot the final result and exact solution
xx = LinRange(0,L,250)

# Create solution trace
exact = scatter(x=xx, y=u_exact(xx,T), mode="lines", name="exact", showlegend=true)

layout = Layout(
    xaxis_title = "x",
    yaxis_title = "u(x,0.5)"
    )

plot([exact, f_euler, b_euler], layout)