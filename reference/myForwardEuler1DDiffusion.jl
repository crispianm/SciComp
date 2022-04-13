# simple forward Euler solver for the 1D heat equation
#   u_t = kappa u_xx  0<x<L, 0<t<T
# with zero-temperature boundary conditions
#   u=0 at x=0,L, t>0
# and prescribed initial temperature
#   u=u_I(x) 0<=x<=L,t=0



using PlotlyJS
# Set problem parameters/functions
kappa = 1.0   # diffusion constant
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

# Plot the final result and exact solution
xx = LinRange(0,L,250)

# Create traces
exact = scatter(x=xx, y=u_exact(xx,T), mode="lines", name="exact", showlegend=true)
num = scatter(x=x, y=u_j, mode="markers", name="num", showlegend=true)

layout = Layout(
    xaxis_title = "x",
    yaxis_title = "u(x,0.5)"
    )

plot([exact, num], layout)