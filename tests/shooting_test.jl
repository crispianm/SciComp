function hopf(u1, u2, beta, sigma=-1.0)

    du1dt = beta*u1 - u2 + sigma*u1*(u1^2 + u2^2)
    du2dt = u1 + beta*u2 + sigma*u2*(u1^2 + u2^2)

    return [du1dt du2dt]
end

function hopf_sol(t, beta, theta=0.0)

    u1 = √(beta) * cos(t + theta)
    u2 = √(beta) * sin(t + theta)

    return [u1 u2]
end