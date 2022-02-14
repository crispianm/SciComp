import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

# function that returns dz/dt
def system(u,t):

    a = 1
    d = 0.1
    b = 0.26

    x = u[0]
    y = u[1]

    dxdt = x*(1-x) - (a*x*y)/(d + x)
    dydt = b*y*(1-(y/x))

    dUdt = np.array([dxdt, dydt])

    return dUdt 

# initial condition
u0 = [1,1]

# time points
t = np.linspace(0,100)

# solve ODE
z = odeint(system,u0,t)

# plot results
plt.plot(t,z[:,0],'b-',label='x')
plt.plot(t,z[:,1],'r--',label='y')
plt.ylabel('response')
plt.xlabel('time')
plt.legend(loc='best')
plt.show()