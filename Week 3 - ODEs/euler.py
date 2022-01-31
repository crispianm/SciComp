import numpy as np
import matplotlib.pyplot as plt

## Functions

def euler_step(f, t, deltat_max, s):
    for i in range(0, len(t) - 1):
        s[i + 1] = s[i] + deltat_max*f(t[i], s[i])

def solve_to():

def solve_ode():
    print('hi')



plt.style.use('seaborn-poster')
%matplotlib inline

# Define parameters
f = lambda t, s: np.exp(-t) # ODE
deltat_max = 0.1 # Step size
t = np.arange(0, 1 + deltat_max, deltat_max) # Numerical grid
s0 = -1 # Initial Condition

# Explicit Euler Method
s = np.zeros(len(t))
s[0] = s0

