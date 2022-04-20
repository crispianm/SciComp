# function algebraic(x, c, arg...)

#     u1, u2 = u
#     du1dt = beta*u1 - u2 + sigma*u1*(u1^2 + u2^2)
#     du2dt = u1 + beta*u2 + sigma*u2*(u1^2 + u2^2)

#     return [du1dt du2dt]
# end


# function hopf2d(u, t, beta, sigma=-1.0, arg...)

#     u1, u2 = u
#     du1dt = beta*u1 - u2 + sigma*u1*(u1^2 + u2^2)
#     du2dt = u1 + beta*u2 + sigma*u2*(u1^2 + u2^2)

#     return [du1dt du2dt]
# end


# function hopf2d_modified(u, t, beta, sigma=-1.0, arg...)

#     u1, u2 = u
#     du1dt = beta*u1 - u2 - sigma*u1*(u1^2 + u2^2) + sigma*u1*(u1^2 + u2^2)^2
#     du2dt = u1 + beta*u2 - sigma*u2*(u1^2 + u2^2) + sigma*u2*(u1^2 + u2^2)^2

#     return [du1dt du2dt]
# end


include("../examples/example_functions.jl")